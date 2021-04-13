#ifndef SM_H
#define SM_H

#include <memory>
#include <chrono>
#include <signal.h>
#include <cmath>
#include <numeric>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <tuple>
#include <functional>

#include "fsm.h"
#include "csm.h"


class SM
{
  public:

    SM();

    SM(
      const unsigned int& max_iterations,
      const unsigned int& num_iterations,
      const unsigned int& start_sample,
      const unsigned int& end_sample,
      const double& xy_uniform_displacement,
      const double& t_uniform_displacement,
      const double& sigma_noise_real,
      const double& sigma_noise_map,
      const double& invalid_rays_randomly_percent,
      const double& invalid_rays_sequentially_percent,
      const unsigned int& size_real_scan,
      const unsigned int& size_map,
      const std::string& method,
      const unsigned int max_counter,
      const unsigned int min_magnification_size,
      const unsigned int max_magnification_size,
      const unsigned int max_recoveries,
      const bool enforce_terminal_constraint,
      const bool enforce_early_gearup);

    ~SM();

  private:

    unsigned int MAX_ITERATIONS;
    unsigned int NUM_ITERATIONS;
    unsigned int START_SAMPLE;
    unsigned int END_SAMPLE;
    double XY_UNIFORM_DISPLACEMENT;
    double T_UNIFORM_DISPLACEMENT;
    double SIGMA_NOISE_REAL;
    double SIGMA_NOISE_MAP;
    double INVALID_RAYS_RANDOMLY_PERCENT;
    double INVALID_RAYS_SEQUENTIALLY_PERCENT;
    unsigned int SIZE_REAL_SCAN;
    unsigned int SIZE_MAP;
    std::string METHOD;
    unsigned int MAX_COUNTER;
    unsigned int MIN_MAGNIFICATION_SIZE;
    unsigned int MAX_MAGNIFICATION_SIZE;
    unsigned int MAX_RECOVERIES;
    bool ENFORCE_TERMINAL_CONSTRAINT;
    bool ENFORCE_EARLY_GEARUP;

    FSMSM::input_params input_params_;

    sm_params input_;
    sm_result output_;

    std::string base_path_;

    std::string real_poses_filename_;
    std::string initial_pose_estimates_;
    std::string pose_errors_filename_;
    std::string position_errors_filename_;
    std::string orientation_errors_filename_;
    std::string trajectories_filename_;
    std::string rotation_times_filename_;
    std::string rotation_iterations_filename_;
    std::string translation_times_filename_;
    std::string translation_iterations_filename_;
    std::string times_filename_;
    std::string intersections_times_filename_;
    std::string num_recoveries_filename_;
    std::string statistics_filename_;
    char dataset_filepath_[];

    // Forward plan (DFT)
    fftw_plan r2rp_;

    // Backward plan (IDFT)
    fftw_plan c2rp_;


    // **** methods

    void cacheFFTW3Plans(const unsigned int& sz);

    std::vector<double> collectErrorBins(
      const std::vector<double>& errors);

    double computePositionError(
      const std::tuple<double,double,double>& pose1,
      const std::tuple<double,double,double>& pose2);
    double computeOrientationError(
      const std::tuple<double,double,double>& pose1,
      const std::tuple<double,double,double>& pose2);
    double computePoseError(
      const std::tuple<double,double,double>& pose1,
      const std::tuple<double,double,double>& pose2);

    void corruptRanges(
      std::vector<double>* scan, const double& sigma);

    void corruptMap(std::vector< std::pair<double,double> >* map,
      const double& sigma);

    void doFinalReport(
      const std::vector< std::vector<double> >& errors_n,
      const std::vector< std::vector<int> >& iterations_n,
      const std::vector< std::vector<double> >& times_n,
      const std::vector< std::vector<double> >& intersections_times_n);

    void doOneReport(
      const std::vector<double>& errors,
      const std::vector<int>& iterations,
      const std::vector<double>& times);

    void initParams(const unsigned int& max_iterations,
      const unsigned int& num_iterations,
      const unsigned int& start_sample,
      const unsigned int& end_sample,
      const double& xy_uniform_displacement,
      const double& t_uniform_displacement,
      const double& sigma_noise_real,
      const double& sigma_noise_map,
      const double& invalid_rays_randomly_percent,
      const double& invalid_rays_sequentially_percent,
      const unsigned int& size_real_scan,
      const unsigned int& size_map,
      const std::string& method,
      const unsigned int max_counter,
      const unsigned int min_magnification_size,
      const unsigned int max_magnification_size,
      const unsigned int max_recoveries,
      const bool enforce_terminal_constraint,
      const bool enforce_early_gearup);

    void initLogs();

    void invalidateRaysRandomly(std::vector<double>* scan,
      const double& percent_invalid);

    void invalidateRaysSequentially(std::vector<double>* scan,
      const double& percent_invalid);

    void generateMapConfiguration(
      const std::tuple<double,double,double>& dataset_pose,
      const std::vector<double>& dataset_scan,
      const bool& do_artificial_360,
      const bool& generate_real_pose,
      std::vector< std::pair<double,double> >* map,
      std::vector<double>* real_scan,
      std::tuple<double,double,double>* real_pose,
      std::vector<double>* virtual_scan,
      std::tuple<double,double,double>* virtual_pose);

    void performTests();

    void testOne(
      const std::vector< double >& real_scan,
      const std::tuple<double,double,double>& real_pose,
      const std::vector< double >& virtual_scan,
      const std::tuple<double,double,double>& virtual_pose,
      const std::vector< std::pair<double,double> >& map,
      const FSMSM::input_params& ip, FSMSM::output_params* op,
      std::tuple<double,double,double>* result_pose);

};

#endif // SM_H
