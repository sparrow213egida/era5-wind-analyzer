#include <iostream>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <numbers>
#include <unordered_map>
#include <vector>
#include <ranges>
#include <fstream>
#include "gdal_priv.h"

struct WindAnalyzer {
    explicit WindAnalyzer(const char *path) {
        /*
         * initializes struct members
         *
         * const char *path - path to the dataset on the user's device
         *
         * if an error occurs while opening the dataset, error is thrown
         * the dataset must contain a short name and valid time for each entry
         */
        GDALAllRegister();

        dataset_unique_ptr_ =
            GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen(path,GA_ReadOnly)));

        if (!dataset_unique_ptr_) {
            throw std::runtime_error("Error occurred while opening the dataset!");
        }

        x_size_ = dataset_unique_ptr_->GetRasterXSize();
        y_size_ = dataset_unique_ptr_->GetRasterYSize();

        dataset_unique_ptr_->GetGeoTransform(dataset_transform_vector_);

        const int raster_cnt = dataset_unique_ptr_->GetRasterCount();
        for (int raster_ind = 1; raster_ind <= raster_cnt; ++raster_ind) {
            GDALRasterBand *raster_band = dataset_unique_ptr_->GetRasterBand(raster_ind);

            const char *short_name = raster_band->GetMetadataItem("GRIB_ELEMENT"); // NOLINT
            if (!short_name) {
                throw std::runtime_error("Error: the dataset does not contain a short name for every entry!");
            }

            const char *valid_time = raster_band->GetMetadataItem("GRIB_VALID_TIME");
            if (!valid_time) {
                throw std::runtime_error("Error: the dataset does not contain a valid time for every entry!");
            }

            if (short_name == std::string("10U")) {
                u_timestamp_to_index_[std::stoll(valid_time) / 3600] = raster_ind;
            } else if (short_name == std::string("10V")) {
                v_timestamp_to_index_[std::stoll(valid_time) / 3600] = raster_ind;
            }
        }
    }

    WindAnalyzer(const WindAnalyzer& other_analyser) = delete;

    WindAnalyzer(WindAnalyzer&& other_analyser) noexcept = default;

    WindAnalyzer& operator=(const WindAnalyzer& other_analyser) = delete;

    WindAnalyzer& operator=(WindAnalyzer&& other_analyser) noexcept = default;

    ~WindAnalyzer() = default;

    [[nodiscard]] std::pair<double, double> GetWindInfo(
        const long long timestamp, const int longitude, const int latitude) {
        /*
         * returns a pair of the speed and the direction of the wind at a particular moment in time
         * in a particular place
         *
         * const long long timestamp - number of the hour, of which the wind info will be calculated (in hours)
         * const int longitude - longitude of the place, where the wind info will be calculated (in degrees)
         * const int latitude - latitude of the place, where the wind info will be calculated (in degrees)
         *
         * return value - pair of speed of the wind (in m/s) and azimuth of the direction of the wind (in degrees)
         *
         * if an error occurs while reading from the dataset, error message is printed and {-1.0, -1.0} is returned
         * if longitude, latitude or timestamp are out of range, an error is also thrown
         */
        if (longitude < -180 || longitude > 180) {
            throw std::out_of_range("Error: longitude must be from -180 to 180!");
        }
        if (latitude < -90 || latitude > 90) {
            throw std::out_of_range("Error: latitude must be from -90 to 90!");
        }

        if (!u_timestamp_to_index_.contains(timestamp) || !v_timestamp_to_index_.contains(timestamp)) {
            throw std::out_of_range("Error: dataset does not contain information, corresponding to this hour!");
        }

        const int u_band_number = u_timestamp_to_index_[timestamp];
        const int v_band_number = v_timestamp_to_index_[timestamp];

        GDALRasterBand *u_band_ptr = dataset_unique_ptr_->GetRasterBand(u_band_number);
        GDALRasterBand *v_band_ptr = dataset_unique_ptr_->GetRasterBand(v_band_number);

        const auto u_ptr = std::make_unique<double>();
        const auto v_ptr = std::make_unique<double>();

        const auto [x_offset, y_offset] = OffsetDegreesToIndexes(longitude, latitude);

        int err_code = u_band_ptr->RasterIO(GF_Read, x_offset, y_offset, 1, 1,
        u_ptr.get(), 1, 1, GDT_Float64, 0, 0);

        if (err_code != 0) {
            std::cerr << "Error occurred while reading from the dataset!" << std::endl;
            return {-1.0, -1.0};
        }

        err_code = v_band_ptr->RasterIO(GF_Read, x_offset, y_offset, 1, 1,
        v_ptr.get(), 1, 1, GDT_Float64, 0, 0);

        if (err_code != 0) {
            std::cerr << "Error occurred while reading from the dataset!" << std::endl;
            return {-1.0, -1.0};
        }

        return {GetWindSpeed(*u_ptr, *v_ptr), GetWindDirection(*u_ptr, *v_ptr)};
    }

protected:
    static double GetWindSpeed(const double u_comp, const double v_comp) {
        /*
         * returns the speed of the wind at a particular moment in time
         *
         * const double u_comp - u-component of the speed of the wind (in m/s)
         * const double v_comp - v-component of the speed of the wind (in m/s)
         *
         * return value - speed of the wind (in m/s)
         */
        return std::sqrt(u_comp * u_comp + v_comp * v_comp);
    }

    static double GetWindDirection(const double u_comp, const double v_comp) {
        /*
         * returns the direction of the wind at a particular moment in time
         *
         * const double u_comp - u-component of the speed of the wind (in m/s)
         * const double v_comp - v-component of the speed of the wind (in m/s)
         *
         * return value - azimuth of the direction of the wind (in degrees)
         */
        const double math_angle = atan2(v_comp, u_comp);

        double azimuth_in_radians;

        if (math_angle > std::numbers::pi / 2) {
            azimuth_in_radians = 5 * std::numbers::pi / 2 - math_angle;
        } else {
            azimuth_in_radians = std::numbers::pi / 2 - math_angle;
        }

        return 180.0 / std::numbers::pi * azimuth_in_radians;
    }

    static std::pair<int, int> SolveForOffsets(const int longitude, const int latitude,
        const GDALGeoTransform& transform_vector) {
        /*
         * solves for x_offset and y_offset using Cramer's method
         *
         * const int longitude - longitude of a place, of which x_offset needs to be calculated
         * const int latitude - latitude of a place, of which y_offset needs to be calculated
         * const GDALGeoTransform& transform_vector - vector of transformation from pixel space to
         * geographic coordinated space
         *
         * return value - pair of x_offset and y_offset for a given point
         */
        const double det = transform_vector[1] * transform_vector[5] - transform_vector[4] * transform_vector[2];
        const double x_det = (longitude - transform_vector[0]) * transform_vector[5] -
            (latitude - transform_vector[3]) * transform_vector[2];
        const double y_det = transform_vector[1] * (latitude - transform_vector[3]) -
            transform_vector[4] * (longitude - transform_vector[0]);

        return {static_cast<int>(std::round(x_det / det)), static_cast<int>(std::round(y_det / det))};
    }

    [[nodiscard]] std::pair<int, int> OffsetDegreesToIndexes(const int longitude, const int latitude) const {
        /*
         * returns a pair of indexes for RasterIO
         *
         * const int longitude - longitude of a given place (in degrees)
         * const int latitude - latitude of a given place (in degrees)
         *
         * return value - a pair of offset indexes, later used for RasterIO
         *
         * if a place with longitude and latitude, passed by the user, has no corresponding entries in
         * the dataset, an error is thrown
         */
        auto [x_offset, y_offset] = SolveForOffsets(longitude, latitude, dataset_transform_vector_);

        if (x_offset < 0 || x_offset > x_size_ || y_offset < 0 || y_offset > y_size_) {
            throw std::out_of_range("Error: the dataset does not contain information about such place!");
        }

        if (x_offset == x_size_) {
            x_offset = x_size_ - 1;
        }
        if (y_offset == y_size_) {
            y_offset = y_size_ - 1;
        }

        return {x_offset, y_offset};
    }

    GDALDatasetUniquePtr dataset_unique_ptr_ = nullptr;
    std::unordered_map<long long, int> u_timestamp_to_index_{};
    std::unordered_map<long long, int> v_timestamp_to_index_{};
    int x_size_ = 0;
    int y_size_ = 0;
    GDALGeoTransform dataset_transform_vector_{};
};

struct WindAnalyzerWithExtraStatistic : WindAnalyzer {
    void WriteMeanWindPowerData(const char *path) {
        /*
         * writes data about windpower in a CSV-file
         *
         * const char *path - path to the CSV-file, in which data will be written
         *
         * no return value
         * if file cannot be opened, an error is thrown
         */
        const std::vector<std::vector<double>> wind_power_data = GetMeanWindPower();

        std::ofstream out(path);
        if (!out.is_open()) {
            std::cerr << "Error: file could not be opened!" << std::endl;
            return;
        }

        out << "Longitude,Latitude,MeanWindpower\n";

        for (int y_ind = 0; y_ind < y_size_; y_ind += 2) {
            for (int x_ind = 0; x_ind < x_size_; x_ind += 2) {
                auto [longitude, latitude] = GetLongitudeAndLatitude(x_ind, y_ind);
                out << longitude << "," << latitude << "," << wind_power_data[y_ind][x_ind] << "\n";
            }
        }

        out.close();
    }
private:
    using WindAnalyzer::WindAnalyzer;

    std::pair<double, double> GetLongitudeAndLatitude(const int x_ind, const int y_ind) const {
        /*
         * returns longitude and latitude from raster indexes
         *
         * const int x_ind - x-coordinate of an index in a raster band
         * const int y_ind - y-coordinate of an index in a raster band
         *
         * return value - pair of longitude and latitude (in degrees)
         */
        const double longitude = dataset_transform_vector_[0] + x_ind * dataset_transform_vector_[1] +
            y_ind * dataset_transform_vector_[2];
        const double latitude = dataset_transform_vector_[3] + x_ind * dataset_transform_vector_[4] +
            y_ind * dataset_transform_vector_[5];

        return {longitude, latitude};
    }

    std::vector<std::vector<double>> GetMeanWindPower() {
        /*
         * returns data about mean windpower
         *
         * has no arguments
         *
         * return value - vector<vector<double>>, where each double value corresponds to
         * mean value of windpower in a corresponding piece of territory
         * if the dataset is not initialized, an error is thrown
         * if the dataset has not enough information for any timestamp, an error is thrown
         */
        if (!dataset_unique_ptr_) {
            std::cerr << "Error: data from the dataset could not be read!" << std::endl;
            return {};
        }

        std::vector<std::vector<double>> mean_power(y_size_);
        for (auto& layer : mean_power) {
            layer.resize(x_size_);
        }

        std::vector<double> cur_u_components(x_size_ * y_size_);
        std::vector<double> cur_v_components(x_size_ * y_size_);

        int cnt_timestamps_taken_into_account = 0;
        for (const long long timestamp : u_timestamp_to_index_ | std::ranges::views::keys) {
            if (v_timestamp_to_index_.contains(timestamp)) {
                UpdateMeanPower(mean_power,
                    u_timestamp_to_index_[timestamp], v_timestamp_to_index_[timestamp],
                    cur_u_components, cur_v_components);
                ++cnt_timestamps_taken_into_account;
            }
        }

        if (cnt_timestamps_taken_into_account == 0) {
            throw std::runtime_error("Error: no matching timestamps found!");
        }

        for (auto& layer : mean_power) {
            for (auto& el : layer) {
                el /= cnt_timestamps_taken_into_account;
            }
        }

        return mean_power;
    }

    void UpdateMeanPower(std::vector<std::vector<double>>& mean_power,
        const int u_band_ind, const int v_band_ind,
        std::vector<double>& u_components, std::vector<double>& v_components) const {
        /*
         * updates value of mean_power (now sum_power) for each raster pixel
         *
         * std::vector<std::vector<double>>& mean power - a vector of sum powers for all raster bands
         * for each raster pixel, which needs to be updated
         * const int u_band_ind - index of raster band, corresponding to u-component of the wind
         * in a given timestamp
         * const int v_band_ind - index of raster band, corresponding to v-component of the wind
         * in a given timestamp
         * std::vector<double>& u_components - vector of u_components for current raster band
         * std::vector<double>& v_components - vector of v_components for current raster band
         *
         * no return value
         *
         * if data from the raster band cannot be read, an error is thrown
         */
        int err_code = dataset_unique_ptr_->GetRasterBand(u_band_ind)->ReadRaster(u_components); // NOLINT
        if (err_code != 0) {
            throw std::runtime_error("Error: could not read from the raster band!");
        }

        err_code = dataset_unique_ptr_->GetRasterBand(v_band_ind)->ReadRaster(v_components); // NOLINT
        if (err_code != 0) {
            throw std::runtime_error("Error: could not read from the raster band!");
        }

        for (int i = 0; i < x_size_ * y_size_; ++i) {
            const double speed = GetWindSpeed(u_components[i], v_components[i]);
            mean_power[i / x_size_][i % x_size_] += speed * speed * speed;
        }
    }
};

int main() {
    constexpr auto path = "/home/sparr/projects/TestTask/full_data.grib";

    auto wind_analyzer = WindAnalyzerWithExtraStatistic(path);
    wind_analyzer.WriteMeanWindPowerData("/home/sparr/projects/TestTask/power_data.csv");
}