#include <sm.h>

/*******************************************************************************
*/
SM::SM(
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
  const bool enforce_early_gearup)
{
  printf("[SM] Starting with custom params \n");
  printf("FSM params: \n");
  printf("max iterations: %u\n", max_iterations);
  printf("num repetitions %u\n", num_iterations);
  printf("(start_sample, end_sample) = (%u,%u)\n", start_sample, end_sample);
  printf("d_xy = +/-%f\n", xy_uniform_displacement);
  printf("d_t = +/-%f\n", t_uniform_displacement);
  printf("sigma_R = %f\n", sigma_noise_real);
  printf("sigma_M = %f\n", sigma_noise_map);
  printf("N_s = %u\n", size_real_scan);
  printf("M_s = %u\n", size_map);
  printf("max_counter = %u\n", max_counter);
  printf("min_magnification_size = %u\n", min_magnification_size);
  printf("max_magnification_size = %u\n", max_magnification_size);
  printf("max_recoveries = %u\n", max_recoveries);
  printf("enforce_terminal_constraint = %d\n", enforce_terminal_constraint);
  printf("enforce_early_gearup = %d\n", enforce_early_gearup);

  // init params
  initParams(
    max_iterations,
    num_iterations,
    start_sample,
    end_sample,
    xy_uniform_displacement,
    t_uniform_displacement,
    sigma_noise_real,
    sigma_noise_map,
    invalid_rays_randomly_percent,
    invalid_rays_sequentially_percent,
    size_real_scan,
    size_map,
    method,
    max_counter,
    min_magnification_size,
    max_magnification_size,
    max_recoveries,
    enforce_terminal_constraint,
    enforce_early_gearup);

  cacheFFTW3Plans(SIZE_REAL_SCAN);


  // init logs
#if defined (LOGS)
  initLogs();
#endif

  // Start tests
  performTests();
}


/*******************************************************************************
*/
SM::~SM()
{
  printf("[SM] Destroying SM\n");
}


/*******************************************************************************
*/
void SM::cacheFFTW3Plans(const unsigned int& sz)
{
  // Create forward  plan
  double* r2r_in;
  double* r2r_out;

  r2r_in = (double*) fftw_malloc(sz * sizeof(double));
  r2r_out = (double*) fftw_malloc(sz * sizeof(double));

  r2rp_ = fftw_plan_r2r_1d(sz, r2r_in, r2r_out, FFTW_R2HC, FFTW_MEASURE);

  // Create backward plan
  fftw_complex* c2r_in;
  double* c2r_out;

  c2r_in = (fftw_complex*) fftw_malloc(sz * sizeof(fftw_complex));
  c2r_out = (double*) fftw_malloc(sz * sizeof(double));

  c2rp_ = fftw_plan_dft_c2r_1d(sz, c2r_in, c2r_out, FFTW_MEASURE);
}


/*******************************************************************************
*/
  std::vector<double>
SM::collectErrorBins(const std::vector<double>& errors)
{
  double b0 = 0.001;
  double b1 = 0.005;
  double b2 = 0.01;
  double b3 = 0.05;

  double num_below_b0 = 0;
  double num_below_b1 = 0;
  double num_below_b2 = 0;
  double num_below_b3 = 0;
  double num_above_b3 = 0;

  for (int i = 0; i < errors.size(); i++)
  {
    if (errors[i] <= b0)
      num_below_b0++;
    else if (errors[i] <= b1)
      num_below_b1++;
    else if (errors[i] <= b2)
      num_below_b2++;
    else if (errors[i] <= b3)
      num_below_b3++;
    else
      num_above_b3++;
  }

  std::vector<double> error_bins;
  error_bins.push_back(num_below_b0);
  error_bins.push_back(num_below_b1);
  error_bins.push_back(num_below_b2);
  error_bins.push_back(num_below_b3);
  error_bins.push_back(num_above_b3);

  return error_bins;
}


/*******************************************************************************
*/
double SM::computePositionError(
  const std::tuple<double,double,double>& pose1,
  const std::tuple<double,double,double>& pose2)
{
  double dx = std::get<0>(pose1) - std::get<0>(pose2);
  double dy = std::get<1>(pose1) - std::get<1>(pose2);

  return sqrt(dx*dx + dy*dy);
}


/*******************************************************************************
*/
double SM::computeOrientationError(
  const std::tuple<double,double,double>& pose1,
  const std::tuple<double,double,double>& pose2)
{
  double dt = std::get<2>(pose1) - std::get<2>(pose2);
  FSMSM::Utils::wrapAngle(&dt);

  return fabs(dt);
}

/*******************************************************************************
*/
double SM::computePoseError(
  const std::tuple<double,double,double>& pose1,
  const std::tuple<double,double,double>& pose2)
{
  double dx = std::get<0>(pose1) - std::get<0>(pose2);
  double dy = std::get<1>(pose1) - std::get<1>(pose2);
  double dt = std::get<2>(pose1) - std::get<2>(pose2);
  FSMSM::Utils::wrapAngle(&dt);

  return sqrt(dx*dx + dy*dy + dt*dt);
}


/*******************************************************************************
*/
  void
SM::corruptRanges(std::vector<double>* scan, const double& sigma)
{
  if (sigma > 0.0)
  {
    std::random_device rand_dev;
    std::mt19937 generator_n(rand_dev());

    std::normal_distribution<double> distribution_n(0.0, sigma);

    for (int i = 0; i < scan->size(); i++)
    {
      double n = 0.0;
      do
      {
        n = distribution_n(generator_n);
      } while (scan->at(i) + n < 0); // ranges should be non-negative

      scan->at(i) += n;
    }
  }
}


/*******************************************************************************
*/
  void
SM::corruptMap(std::vector< std::pair<double,double> >* map,
  const double& sigma)
{
  if (sigma > 0.0)
  {
    std::random_device rand_dev_x;
    std::random_device rand_dev_y;
    std::mt19937 generator_x(rand_dev_x());
    std::mt19937 generator_y(rand_dev_y());

    std::normal_distribution<double> distribution_x(0.0, sigma);
    std::normal_distribution<double> distribution_y(0.0, sigma);

    for (int i = 0; i < map->size(); i++)
    {
      double dx = distribution_x(generator_x);
      map->at(i).first += dx;
      double dy = distribution_y(generator_y);
      map->at(i).second += dy;
    }
  }
}


/*******************************************************************************
*/
  void
SM::doFinalReport(
  const std::vector< std::vector<double> >& errors_n,
  const std::vector< std::vector<int> >& iterations_n,
  const std::vector< std::vector<double> >& times_n,
  const std::vector< std::vector<double> >& intersections_times_n)
{
  printf("________________________________________________________________________\n");
  printf("Displacement uniform in (+/-%.2f,+/-%.2f,+/-%.3f)\n",
    XY_UNIFORM_DISPLACEMENT, XY_UNIFORM_DISPLACEMENT, T_UNIFORM_DISPLACEMENT);
  printf("________________________________________________________________________\n");
  printf("Noise acting on laser's range measurements: Normal with sigma = %f\n",
    SIGMA_NOISE_REAL);
  printf("Noise acting on map's range measurements: Normal with sigma = %f\n",
    SIGMA_NOISE_MAP);
  printf("Invalidated rays in random percentage = %f\n", INVALID_RAYS_RANDOMLY_PERCENT);
  printf("Invalidated rays in sequence percentage = %f\n", INVALID_RAYS_SEQUENTIALLY_PERCENT);
  printf("________________________________________________________________________\n");

  double b0 = 0.001;
  double b1 = 0.005;
  double b2 = 0.01;
  double b3 = 0.05;

  double num_below_b0 = 0;
  double num_below_b1 = 0;
  double num_below_b2 = 0;
  double num_below_b3 = 0;
  double num_above_b3 = 0;

  double mean_error = 0.0;
  double mean_iterations = 0.0;
  double mean_time = 0.0;
  double mean_intersections_time = 0.0;

  for (int i = 0; i < errors_n.size(); i++)
  {
    std::vector<double> error_bins = collectErrorBins(errors_n[i]);
    num_below_b0 += error_bins[0];
    num_below_b1 += error_bins[1];
    num_below_b2 += error_bins[2];
    num_below_b3 += error_bins[3];
    num_above_b3 += error_bins[4];

    for(int k = 0; k < errors_n[i].size(); k++)
    {
      mean_error += errors_n[i][k];
      mean_iterations += iterations_n[i][k];
      mean_time += times_n[i][k];
      mean_intersections_time += intersections_times_n[i][k];
    }
  }

  size_t denom = (errors_n.size() * errors_n[0].size());

  printf("Below %.3f: %f (%.0f/%zu)\n",
    b0, 100*num_below_b0 / denom, num_below_b0, denom);
  printf("Below %.3f: %f (%.0f/%zu)\n",
    b1, 100*num_below_b1 / denom, num_below_b1, denom);
  printf("Below %.3f: %f (%.0f/%zu)\n",
    b2, 100*num_below_b2 / denom, num_below_b2, denom);
  printf("Below %.3f: %f (%.0f/%zu)\n",
    b3, 100*num_below_b3 / denom, num_below_b3, denom);
  printf("Above %.3f: %f (%.0f/%zu)\n",
    b3, 100*num_above_b3 / denom, num_above_b3, denom);

  printf("_________________________________\n");
  printf("Final mean error: %f\n", mean_error / denom);
  printf("Final mean number of iterations: %f\n", mean_iterations / denom);
  printf("Final mean execution time: %f\n", mean_time / denom);
  printf("Final mean intersection-finding time: %f\n", mean_intersections_time / denom);
  printf("_________________________________\n");

#if defined (LOGS)
  std::ofstream statistics_file(statistics_filename_.c_str(), std::ios::app);
  if (statistics_file.is_open())
  {
    statistics_file << "Displacement uniform in ["
      << "+/-" << XY_UNIFORM_DISPLACEMENT << ","
      << "+/-" << XY_UNIFORM_DISPLACEMENT << ","
      << "[" << -T_UNIFORM_DISPLACEMENT << ","
      << T_UNIFORM_DISPLACEMENT <<
      "]]"<< std::endl;
    statistics_file << "Noise acting on laser's range measurements: Normal with sigma = "
      << SIGMA_NOISE_REAL << std::endl;
    statistics_file << "Noise acting on map's range measurements: Normal with sigma = "
      << SIGMA_NOISE_MAP << std::endl;
    statistics_file << "Invalidated rays randomly percentage = "
      << INVALID_RAYS_RANDOMLY_PERCENT << std::endl;
    statistics_file << "Invalidated rays sequentially percentage = "
      << INVALID_RAYS_SEQUENTIALLY_PERCENT << std::endl;
    statistics_file << "Below " << b0 << ": " << 100*num_below_b0 / denom
      << "(" << num_below_b0 << "/" << denom << ")" << std::endl;
    statistics_file << "Below " << b1 << ": " << 100*num_below_b1 / denom
      << "(" << num_below_b1 << "/" << denom << ")" << std::endl;
    statistics_file << "Below " << b2 << ": " << 100*num_below_b2 / denom
      << "(" << num_below_b2 << "/" << denom << ")" << std::endl;
    statistics_file << "Below " << b3 << ": " << 100*num_below_b3 / denom
      << "(" << num_below_b3 << "/" << denom << ")" << std::endl;
    statistics_file << "Above " << b3 << ": " << 100*num_above_b3 / denom
      << "(" << num_above_b3 << "/" << denom << ")" << std::endl;
    statistics_file << "Mean error: "
      << mean_error / denom << " m" << std::endl;
    statistics_file << "Mean number of iterations: "
      << mean_iterations / denom << std::endl;
    statistics_file << "Mean execution time: "
      << mean_time / denom << " sec" << std::endl;
    statistics_file.close();
  }
#endif
}


/*******************************************************************************
*/
  void
SM::doOneReport(
  const std::vector<double>& errors,
  const std::vector<int>& iterations,
  const std::vector<double>& times)
{
  printf("__________________________________________________\n");
  printf("Displacement uniform in (+/-%.2f,+/-%.2f,+/-%.3f)\n",
    XY_UNIFORM_DISPLACEMENT, XY_UNIFORM_DISPLACEMENT, T_UNIFORM_DISPLACEMENT);
  printf("_________________________________\n");

  std::vector<double> error_bins = collectErrorBins(errors);

  size_t denom = (errors.size());

  double b0 = 0.001;
  double b1 = 0.005;
  double b2 = 0.01;
  double b3 = 0.05;

  printf("Below %.3f: %f (%.0f/%zu)\n",
    b0, 100*error_bins[0] / denom, error_bins[0], denom);

  if (error_bins[1] > 0)
  {
    printf("\033[0;31m");
    printf("Below %.3f: %f (%.0f/%zu)\n",
      b1, 100*error_bins[1] / denom, error_bins[1], denom);
    printf("\033[0m");
  }
  else
    printf("Below %.3f: %f (%.0f/%zu)\n",
      b1, 100*error_bins[1] / denom, error_bins[1], denom);

  if (error_bins[2] > 0)
  {
    printf("\033[0;31m");
    printf("Below %.3f: %f (%.0f/%zu)\n",
      b2, 100*error_bins[2] / denom, error_bins[2], denom);
    printf("\033[0m");
  }
  else
    printf("Below %.3f: %f (%.0f/%zu)\n",
      b2, 100*error_bins[2] / denom, error_bins[2], denom);

  if (error_bins[3] > 0)
  {
    printf("\033[0;31m");
    printf("Below %.3f: %f (%.0f/%zu)\n",
      b3, 100*error_bins[3] / denom, error_bins[3], denom);
    printf("\033[0m");
  }
  else
    printf("Below %.3f: %f (%.0f/%zu)\n",
      b3, 100*error_bins[3] / denom, error_bins[3], denom);

  if (error_bins[4] > 0)
  {
    printf("\033[0;31m");
    printf("Above %.3f: %f (%.0f/%zu)\n",
      b3, 100*error_bins[4] / denom, error_bins[4], denom);
    printf("\033[0m");
  }
  else
    printf("Above %.3f: %f (%.0f/%zu)\n",
      b3, 100*error_bins[4] / denom, error_bins[4], denom);

  printf("_________________________________\n");
  printf("Mean error: %f\n",
    accumulate(errors.begin(), errors.end(), 0.0) / errors.size());
  printf("Mean number of iterations: %f\n",
    static_cast<double>(accumulate(iterations.begin(), iterations.end(), 0.0))
    / iterations.size());
  printf("Mean execution time: %f\n",
    accumulate(times.begin(), times.end(), 0.0) / times.size());
}


/*******************************************************************************
*/
  void
SM::generateMapConfiguration(
  const std::tuple<double,double,double>& dataset_pose,
  const std::vector<double>& dataset_scan,
  const bool& do_artificial_360,
  const bool& generate_real_pose,
  std::vector< std::pair<double,double> >* map,
  std::vector<double>* real_scan,
  std::tuple<double,double,double>* real_pose,
  std::vector<double>* virtual_scan,
  std::tuple<double,double,double>* virtual_pose)
{
  map->clear();
  real_scan->clear();

  // CONSTRUCT MAP -------------------------------------------------------------
  // How to complete the scan if it doesn't range over 2π?
  if (do_artificial_360)
  {
    std::tuple<double,double,double> origin_of_generated_pose;
    std::tuple<double,double,double> map_origin;

    int completion_fashion = 4;
    if (completion_fashion != 5)
    {
      *real_scan = dataset_scan;
      FSMSM::ScanCompletion::completeScan(real_scan, completion_fashion);
      FSMSM::Utils::scan2points(*real_scan, dataset_pose, map);
      *real_pose = dataset_pose;
      origin_of_generated_pose = dataset_pose;
    }
    else
    {
      // The real scan in points in 2D, aka the map
      FSMSM::ScanCompletion::completeScan5(dataset_pose, dataset_scan,
        pow(2,0)*dataset_scan.size(), real_scan, map, &map_origin);
      *real_pose = map_origin;
      origin_of_generated_pose = map_origin;
    }

    *map = FSMSM::X::findExact(origin_of_generated_pose, *map, SIZE_MAP);
  }
  else // just impress the rays onto the plane over 2π
  {
    // 361-1 = 360
    std::vector<double> map_scan = dataset_scan;
    map_scan.erase(map_scan.begin() + map_scan.size()-1);

    FSMSM::Utils::scan2points(map_scan, dataset_pose, map);

    // This is the map
    *map = FSMSM::X::findExact(dataset_pose, *map, SIZE_MAP);
  }


  // GENERATE REAL POSE (perhaps, sir?) ----------------------------------------
  // Generate a new pose in the map or take that from the dataset
  if (generate_real_pose)
  {
    //while(!Utils::generatePose(dataset_pose, *map,
    //10*XY_UNIFORM_DISPLACEMENT, 0.0, 0.25,
    //real_pose));

    while(!FSMSM::Utils::generatePoseWithinMap(*map, 0.5, real_pose));

    // Limit precision
    // http://www.cplusplus.com/forum/general/222965/#msg1022316
    double f = pow(10, 6);
    std::get<0>(*real_pose) = ((int)(std::get<0>(*real_pose)*f))/f;
    std::get<1>(*real_pose) = ((int)(std::get<1>(*real_pose)*f))/f;
    std::get<2>(*real_pose) = ((int)(std::get<2>(*real_pose)*f))/f;

    // sample 96 debugging
    //std::get<0>(*real_pose) = 7.101889;
    //std::get<1>(*real_pose) = 2.048786;
    //std::get<2>(*real_pose) = 2.408916;
  }
  else
    *real_pose = dataset_pose;

  // Set orientation to 0.0 for quick visual inspection of angular errors.
  // It doesn't matter anyway, this is an independent variable
  // TODO comment-out
  //std::get<2>(*real_pose) = 0.0;std::get<2>(dataset_pose);

  // GENERATE VIRTUAL POSE -----------------------------------------------------
  while(!FSMSM::Utils::generatePose(*real_pose, *map,
      XY_UNIFORM_DISPLACEMENT, T_UNIFORM_DISPLACEMENT, 0.25, virtual_pose));

  std::vector< std::pair<double,double> > v_intersections =
    FSMSM::X::findExact(*virtual_pose, *map, SIZE_REAL_SCAN);
  FSMSM::Utils::points2scan(v_intersections, *virtual_pose, virtual_scan);

  //std::get<0>(*virtual_pose) = 7.174819;
  //std::get<1>(*virtual_pose) = 2.049089;
  //std::get<2>(*virtual_pose) = 2.890520;

  // New real pose means new real scan
  // COMPUTE REAL SCAN ---------------------------------------------------------
  std::vector< std::pair<double,double> > intersections =
    FSMSM::X::findExact(*real_pose, *map, SIZE_REAL_SCAN);
  FSMSM::Utils::points2scan(intersections, *real_pose, real_scan);

#if defined (STORE)
  std::vector< std::pair<double,double> > v_intersections =
    FSMSM::X::findExact(*virtual_pose, *map, real_scan->size());
  std::vector<double> virtual_scan;
  FSMSM::Utils::points2scan(v_intersections, *virtual_pose, &virtual_scan);

  Dump::map(*map, base_path_ + "/../matlab/map_dump.m");
  Dump::scan(*real_scan, *real_pose, virtual_scan, *virtual_pose,
    base_path_ + "/../matlab/scan_dump.m");
#endif
}




/*******************************************************************************
*/
  void
SM::initLogs()
{
  std::string configuration_str =
    std::to_string(NUM_ITERATIONS) + "_"
    + std::to_string(MAX_ITERATIONS) + "_"
    + std::to_string(XY_UNIFORM_DISPLACEMENT) + "_"
    + std::to_string(T_UNIFORM_DISPLACEMENT) + "_"
    + std::to_string(SIGMA_NOISE_REAL) + "_"
    + std::to_string(SIGMA_NOISE_MAP) + "_"
    + std::to_string(INVALID_RAYS_RANDOMLY_PERCENT) + "_"
    + std::to_string(INVALID_RAYS_SEQUENTIALLY_PERCENT);

  std::string base_log_path_str =
    base_path_ + "/../logs/" + METHOD + "/" + configuration_str;

  // Log real poses to file ----------------------------------------------------
  real_poses_filename_ =
    base_log_path_str + "/poses/real_poses.txt";

  // Log virtual starting poses to file ----------------------------------------
  initial_pose_estimates_ =
    base_log_path_str + "/poses/initial_pose_estimates.txt";

  // Log pose errors to file ---------------------------------------------------
  pose_errors_filename_ =
    base_log_path_str + "/errors/pose_errors.txt";

  // Log position errors to file -----------------------------------------------
  position_errors_filename_ =
    base_log_path_str + "/errors/position_errors.txt";

  // Log orientation errors to file --------------------------------------------
  orientation_errors_filename_ =
    base_log_path_str + "/errors/orientation_errors.txt";

  // Log trajectories to file --------------------------------------------------
  trajectories_filename_ =
    base_log_path_str + "/trajectories/trajectories.txt";

  // Log rotation execution times to file --------------------------------------
  rotation_times_filename_ =
    base_log_path_str + "/times/rotation_times.txt";

  // Log translation execution times to file -----------------------------------
  translation_times_filename_ =
    base_log_path_str + "/times/translation_times.txt";

  // Log intersections times to file -------------------------------------------
  intersections_times_filename_ =
    base_log_path_str + "/times/intersections_times.txt";

  // Log total execution times to file -----------------------------------------
  times_filename_ =
    base_log_path_str + "/times/total_times.txt";

  // Log rotation iterations to file -------------------------------------------
  rotation_iterations_filename_ =
    base_log_path_str + "/iterations/rotation_iterations.txt";

  // Log translation iterations to file ----------------------------------------
  translation_iterations_filename_ =
    base_log_path_str + "/iterations/translation_iterations.txt";

  // Log number of recoveries to file ------------------------------------------
  num_recoveries_filename_ =
    base_log_path_str + "/recoveries/num_recoveries.txt";

  // Log statistics to file ----------------------------------------------------
  statistics_filename_ =
    base_log_path_str + "/statistics/statistics.txt";

  std::ofstream real_poses_file(real_poses_filename_.c_str());
  std::ofstream virtual_starting_poses_file(initial_pose_estimates_.c_str());
  std::ofstream pose_errors_file(pose_errors_filename_.c_str());
  std::ofstream position_errors_file(position_errors_filename_.c_str());
  std::ofstream orientation_errors_file(orientation_errors_filename_.c_str());
  std::ofstream trajectories_file(trajectories_filename_.c_str());
  std::ofstream rotation_times_file(rotation_times_filename_.c_str());
  std::ofstream translation_times_file(translation_times_filename_.c_str());
  std::ofstream intersections_times_file(intersections_times_filename_.c_str());
  std::ofstream times_file(times_filename_.c_str());
  std::ofstream rotation_iterations_file(rotation_iterations_filename_.c_str());
  std::ofstream translation_iterations_file(translation_iterations_filename_.c_str());
  std::ofstream num_recoveries_file(num_recoveries_filename_.c_str());
  std::ofstream statistics_file(statistics_filename_.c_str());
}


/*******************************************************************************
*/
  void
SM::initParams(const unsigned int& max_iterations,
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
  const bool enforce_early_gearup)
{
  MAX_ITERATIONS = max_iterations;
  NUM_ITERATIONS = num_iterations;
  START_SAMPLE = start_sample;
  END_SAMPLE = end_sample;
  XY_UNIFORM_DISPLACEMENT  = xy_uniform_displacement;
  T_UNIFORM_DISPLACEMENT = t_uniform_displacement;
  SIGMA_NOISE_REAL = sigma_noise_real;
  SIGMA_NOISE_MAP = sigma_noise_map;
  INVALID_RAYS_RANDOMLY_PERCENT = invalid_rays_randomly_percent;
  INVALID_RAYS_SEQUENTIALLY_PERCENT = invalid_rays_sequentially_percent;
  SIZE_REAL_SCAN = size_real_scan;
  SIZE_MAP = size_map;
  METHOD = method;
  MAX_COUNTER = max_counter;
  MIN_MAGNIFICATION_SIZE = min_magnification_size;
  MAX_MAGNIFICATION_SIZE = max_magnification_size;
  MAX_RECOVERIES = max_recoveries;
  ENFORCE_TERMINAL_CONSTRAINT = enforce_terminal_constraint;
  ENFORCE_EARLY_GEARUP = enforce_early_gearup;

  input_params_.num_iterations = max_iterations;
  input_params_.xy_bound = xy_uniform_displacement;
  input_params_.t_bound = t_uniform_displacement;
  input_params_.sigma_noise_real = SIGMA_NOISE_REAL;
  input_params_.sigma_noise_map = SIGMA_NOISE_MAP;
  input_params_.max_counter = max_counter;
  input_params_.min_magnification_size = min_magnification_size;
  input_params_.max_magnification_size = max_magnification_size;
  input_params_.max_recoveries = max_recoveries;
  input_params_.enforce_terminal_constraint = enforce_terminal_constraint;
  input_params_.enforce_early_gearup = enforce_early_gearup;

  std::string cwd("\0",FILENAME_MAX+1);
  base_path_ = getcwd(&cwd[0],cwd.capacity());
}


/*******************************************************************************
*/
  void
SM::invalidateRaysRandomly(std::vector<double>* scan,
  const double& percent_invalid)
{
  if (percent_invalid > 0.0)
  {
    std::random_device rand_dev;
    std::mt19937 generator_i(rand_dev());
    std::uniform_real_distribution<double> distribution_i(0,1);

    for (int i = 0; i < scan->size(); i++)
    {
      double ii = distribution_i(generator_i);
      if (ii <= percent_invalid)
        scan->at(i) = 0.0;
    }
  }
}


/*******************************************************************************
*/
  void
SM::invalidateRaysSequentially(std::vector<double>* scan,
  const double& percent_invalid)
{
  if (percent_invalid > 0.0)
  {
    // The number of invalid rays
    int num_invalid = floor(percent_invalid * scan->size());

    // Throw a dice to see from which ray forward rays are going to be invalidated
    std::random_device rand_dev;
    std::mt19937 generator_i(rand_dev());
    std::uniform_real_distribution<double> distribution_i(0, scan->size()-1);

    int start = floor(distribution_i(generator_i));
    int end = 0;

    if (start + num_invalid < scan->size())
      end = start + num_invalid;
    else
      end = scan->size();

    for (int i = start; i < end; i++)
      scan->at(i) = 0.0;
  }
}


/*******************************************************************************
*/
void SM::performTests()
{
  // The base string of one dataset
  std::string base_str = base_path_ + "/../dataset/dataset_";

  // The errors, iterations, and execution times per dataset,
  // for all iterations
  std::vector< std::vector<double> > errors_n;
  std::vector< std::vector<int> > iterations_n;
  std::vector< std::vector<double> > times_n;
  std::vector< std::vector<double> > intersections_times_n;

  for (int n = 0; n < NUM_ITERATIONS; n++)
  {
    printf("*                               *\n");
    printf("*********************************\n");
    printf("*                               *\n");
    printf("*                               *\n");
    printf("*                               *\n");
    printf("          ITERATION #%d          \n", n);
    printf("*                               *\n");
    printf("*                               *\n");
    printf("*                               *\n");
    printf("*********************************\n");
    printf("*                               *\n");

    // The errors, iterations, and execution times per dataset, per iteration
    std::vector<double> pose_errors;
    std::vector<int> iterations;
    std::vector<double> times;
    std::vector<double> intersections_times;

    for (int i = START_SAMPLE; i < END_SAMPLE; i++)
    {
      printf("_________________________________\n");
      printf("Processing sample no. #%d\n", i);

      // The resulting pose, number of execution iterations, and execution time
      // for the i-th sample of the dataset
      std::tuple<double,double,double> result_pose;

      FSMSM::output_params output_params_;

      // Read dataset i
      std::string dataset_filepath_str = base_str + std::to_string(i) + ".txt";
      int l = dataset_filepath_str.length();
      char dataset_filepath_char[l+1];
      strcpy(dataset_filepath_char, dataset_filepath_str.c_str());

      // Read scan and the pose from which it was taken from the dataset
      std::vector<double> dataset_real_scan;
      std::tuple<double,double,double> dataset_real_pose;
      FSMSM::DatasetUtils::dataset2rangesAndPose(dataset_filepath_char,
        &dataset_real_scan, &dataset_real_pose);

      // Get the map, the real scan, the real pose, and the virtual pose
      std::vector< std::pair<double,double> > true_map;
      std::vector<double> real_scan;
      std::vector<double> virtual_scan;
      std::tuple<double,double,double> real_pose;
      std::tuple<double,double,double> virtual_pose;
      bool do_artificial_360 = true;
      bool generate_real_pose = true;

      generateMapConfiguration(dataset_real_pose, dataset_real_scan,
        do_artificial_360, generate_real_pose,
        &true_map, &real_scan, &real_pose, &virtual_scan, &virtual_pose);

      // For design/debugging purposes
#if defined(TEST_ROTATION_ONLY_CONT) || defined(TEST_ROTATION_ONLY_DISC)
      std::get<0>(virtual_pose) = std::get<0>(real_pose);
      std::get<1>(virtual_pose) = std::get<1>(real_pose);

#if defined (TEST_ROTATION_ONLY_DISC)
      double ang_inc = 2*M_PI/real_scan.size();
      std::get<2>(virtual_pose) = std::get<2>(real_pose) + 49*ang_inc;
#endif
#endif

#if defined (TEST_TRANSLATION_ONLY)
      std::get<2>(virtual_pose) = std::get<2>(real_pose);
#endif

      // Dump initial conf
#if defined (STORE)
      std::vector< std::pair<double,double> > r;
      r = FSMSM::X::findExact(real_pose, true_map, real_scan.size());
      std::vector< std::pair<double,double> > v;
      v = FSMSM::X::findExact(virtual_pose, true_map, real_scan.size());
      Dump::points(r,v,0, base_path_ + "/../matlab/points_dump.m");
      Dump::map(true_map, base_path_ + "/../matlab/map_dump.m");
#endif

      // -----------------------------------------------------------------------
      // This is the meat ------------------------------------------------------
      testOne(real_scan, real_pose, virtual_scan, virtual_pose, true_map,
        input_params_, &output_params_, &result_pose);
      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------

      // Calculate pose error
      double pose_error = computePoseError(real_pose, result_pose);
      double position_error = computePositionError(real_pose, result_pose);
      double orientation_error = computeOrientationError(real_pose, result_pose);

      if (pose_error > 0.001)
      {
        printf("\033[0;31m");
        printf("SSAMPLE NO. %d\n", i);
        printf("\033[0m");
      }

#if defined (DEBUG)
      printf("Position error = %f\n", position_error);
      printf("Orientation error = %f\n", orientation_error);
      printf("Pose error = %f\n", pose_error);
#endif

      // Log results for i-th dataset
      pose_errors.push_back(pose_error);
      iterations.push_back(output_params_.translation_iterations);
      times.push_back(output_params_.exec_time);
      intersections_times.push_back(output_params_.intersections_times);

#if defined (LOGS)
      // Log real pose
      std::ofstream real_poses_file(
        real_poses_filename_.c_str(), std::ios::app);
      if (real_poses_file.is_open())
      {
        real_poses_file <<
          std::get<0>(real_pose) << ", " <<
          std::get<1>(real_pose) << ", " <<
          std::get<2>(real_pose) << std::endl;
        real_poses_file.close();
      }

      // Log initial pose estimate
      std::ofstream initial_pose_estimates_file(
        initial_pose_estimates_.c_str(), std::ios::app);
      if (initial_pose_estimates_file.is_open())
      {
        initial_pose_estimates_file <<
          std::get<0>(virtual_pose) << ", " <<
          std::get<1>(virtual_pose) << ", " <<
          std::get<2>(virtual_pose) << std::endl;
        initial_pose_estimates_file.close();
      }

      // Log pose error
      std::ofstream pose_errors_file(
        pose_errors_filename_.c_str(), std::ios::app);
      if (pose_errors_file.is_open())
      {
        pose_errors_file << pose_error << std::endl;
        pose_errors_file.close();
      }

      // Log position error
      std::ofstream position_errors_file(
        position_errors_filename_.c_str(), std::ios::app);
      if (position_errors_file.is_open())
      {
        position_errors_file << position_error << std::endl;
        position_errors_file.close();
      }

      // Log orientation error
      std::ofstream orientation_errors_file(
        orientation_errors_filename_.c_str(), std::ios::app);
      if (orientation_errors_file.is_open())
      {
        orientation_errors_file << orientation_error << std::endl;
        orientation_errors_file.close();
      }

      // Log trajectories
      std::ofstream trajectories_file(
        trajectories_filename_.c_str(), std::ios::app);
      if (trajectories_file.is_open())
      {
        for (unsigned int i = 0; i < output_params_.trajectory.size(); i++)
        {
          trajectories_file << std::get<0>(output_params_.trajectory[i]) << ", "
            << std::get<1>(output_params_.trajectory[i]) << ", "
            << std::get<2>(output_params_.trajectory[i]) << std::endl;
        }
        trajectories_file.close();
      }

      // Log rotation time
      std::ofstream rotation_times_file(
        rotation_times_filename_.c_str(), std::ios::app);
      if (rotation_times_file.is_open())
      {
        rotation_times_file << output_params_.rotation_times << std::endl;
        rotation_times_file.close();
      }

      // Log translation time
      std::ofstream translation_times_file(
        translation_times_filename_.c_str(), std::ios::app);
      if (translation_times_file.is_open())
      {
        translation_times_file << output_params_.translation_times << std::endl;
        translation_times_file.close();
      }

      // Log intersections time
      std::ofstream intersections_times_file(
        intersections_times_filename_.c_str(), std::ios::app);
      if (intersections_times_file.is_open())
      {
        intersections_times_file << output_params_.intersections_times << std::endl;
        intersections_times_file.close();
      }

      // Log rotation iterations
      std::ofstream rotation_iterations_file(
        rotation_iterations_filename_.c_str(), std::ios::app);
      if (rotation_iterations_file.is_open())
      {
        rotation_iterations_file <<
          output_params_.rotation_iterations << std::endl;
        rotation_iterations_file.close();
      }

      // Log translation iterations
      std::ofstream translation_iterations_file(
        translation_iterations_filename_.c_str(), std::ios::app);
      if (translation_iterations_file.is_open())
      {
        translation_iterations_file <<
          output_params_.translation_iterations << std::endl;
        translation_iterations_file.close();
      }

      // Log number of recoveries
      std::ofstream num_recoveries_file(
        num_recoveries_filename_.c_str(), std::ios::app);
      if (num_recoveries_file.is_open())
      {
        num_recoveries_file << output_params_.num_recoveries << std::endl;
        num_recoveries_file.close();
      }

      // Log total execution time
      std::ofstream times_file(times_filename_.c_str(), std::ios::app);
      if (times_file.is_open())
      {
        times_file << output_params_.exec_time << std::endl;
        times_file.close();
      }
#endif

#if defined (REPORTS)
      doOneReport(pose_errors, iterations, times);
#endif
    }

    errors_n.push_back(pose_errors);
    iterations_n.push_back(iterations);
    times_n.push_back(times);
    intersections_times_n.push_back(intersections_times);

    printf("************************************************************************\n");
    printf("                   Final report over %d/%d iterations                   \n",
      n+1, NUM_ITERATIONS);
    printf("************************************************************************\n");
    doFinalReport(errors_n, iterations_n, times_n, intersections_times_n);
    printf("************************************************************************\n");
  }
}


/*******************************************************************************
 * Corrupt the real scan first, the map second. Test third
 */
  void
SM::testOne(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& real_pose,
  const std::vector< double >& virtual_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const FSMSM::input_params& ip, FSMSM::output_params* op,
  std::tuple<double,double,double>* result_pose)
{
  // ---------------------------------------------------------------------------
  // real scan stuff
  // Are the real scan's ranges are corrupted by noise?
  // ---------------------------------------------------------------------------
  std::vector<double> real_scan_corrupted = real_scan;
  corruptRanges(&real_scan_corrupted, SIGMA_NOISE_REAL);

  // Invalidate some rays
  invalidateRaysRandomly(&real_scan_corrupted, INVALID_RAYS_RANDOMLY_PERCENT);

  // Invalidate a whole sequential region
  invalidateRaysSequentially(&real_scan_corrupted,
    INVALID_RAYS_SEQUENTIALLY_PERCENT);

  // ---------------------------------------------------------------------------
  // virtual scan stuff
  // Are the virtual scan's ranges are corrupted by noise?
  // ---------------------------------------------------------------------------
  std::vector<double> virtual_scan_corrupted = virtual_scan;
  corruptRanges(&virtual_scan_corrupted, SIGMA_NOISE_REAL);

  std::vector< std::pair<double,double> > virtual_scan_points;
  FSMSM::Utils::scan2points(virtual_scan_corrupted, virtual_pose,
    &virtual_scan_points);

  // ---------------------------------------------------------------------------
  // map stuff
  // Is the map corrupted by noise?
  // ---------------------------------------------------------------------------
  std::vector< std::pair<double,double> > map_corrupted = map;
  corruptMap(&map_corrupted, SIGMA_NOISE_MAP);

  // SCAN-MATCHING: instead of the map you get the virtual scan points

#if defined (PRINTS)
  printf("real pose   (%f,%f,%f)\n",
    std::get<0>(real_pose),
    std::get<1>(real_pose),
    std::get<2>(real_pose));
#endif

  // Test per method
  if (METHOD.compare("FSM") == 0)
  {
    FSMSM::Match::fmtdbh(real_scan_corrupted, virtual_pose, virtual_scan_points,
      METHOD, r2rp_, c2rp_, ip, op, result_pose);
  }

#if defined (PRINTS)
  printf("______________________________________________________________\n");
  printf("input pose  (%f,%f,%f)\n",
    std::get<0>(virtual_pose),
    std::get<1>(virtual_pose),
    std::get<2>(virtual_pose));
  printf("real pose   (%f,%f,%f)\n",
    std::get<0>(real_pose),
    std::get<1>(real_pose),
    std::get<2>(real_pose));
  printf("output pose (%f,%f,%f)\n",
    std::get<0>(*result_pose),
    std::get<1>(*result_pose),
    std::get<2>(*result_pose));
  printf("______________________________________________________________\n");
#endif
}
