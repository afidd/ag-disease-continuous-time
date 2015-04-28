#include <string>
#include <sstream>
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "disease.hpp"
#include "scenario.hpp"
#include "hdf_file.hpp"
#include "ensemble.hpp"
#include "naadsm_xml.hpp"
#include "contact_version.hpp"




/*! Save the whole trajectory.
 */
class TrajectorySave : public TrajectoryObserver
{
  std::vector<TrajectoryEntry> trajectory_;
 public:
  virtual void Step(TrajectoryEntry sirt) override {
    trajectory_.emplace_back(sirt);
  }
  virtual const std::vector<TrajectoryEntry>& Trajectory() const {
    return trajectory_; }
};



/*! Save the trajectory every time any of SIR change by a percent.
 */
class PercentTrajectorySave : public TrajectoryObserver
{
  int64_t step_{0};
  int64_t threshhold_{0};
  double percent_{0.0001};

  TrajectoryEntry last_{0,0,0,0};
  std::vector<TrajectoryEntry> trajectory_;
 public:
  PercentTrajectorySave() {}

  virtual void Step(TrajectoryEntry sirt) override {
    if (0==step_) {
      last_=sirt;
      threshhold_=std::floor(percent_*(sirt.s+sirt.i+sirt.r));
      trajectory_.emplace_back(sirt);
    } else {
      bool ps=std::abs(sirt.s-last_.s)>threshhold_;
      bool pi=std::abs(sirt.i-last_.i)>threshhold_;
      bool pr=std::abs(sirt.r-last_.r)>threshhold_;
      if (ps||pi||pr) {
        trajectory_.emplace_back(sirt);
        last_=sirt;
      }
    }
    ++step_;
  }
  virtual const std::vector<TrajectoryEntry>& Trajectory() const {
    return trajectory_; }
};




int main(int argc, char *argv[]) {
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR with demographics.");
  int64_t individual_cnt=100;
  int64_t infected_cnt=1;
  int64_t recovered_cnt=0;

  int run_cnt=1;
  size_t rand_seed=1;
  // Time is in years.
  std::vector<Parameter> parameters;
  parameters.emplace_back(Parameter{ADParam::FirstFarm, "firstfarm", 0,
    "first farm infected"});
  double end_time=365.0;
  int thread_cnt=1;
  std::string log_level;
  std::string data_file("run.h5");
  bool save_file=false;
  std::string translation_file;
  std::string scenario_file;
  std::string herd_file;
  bool test=true;

  desc.add_options()
    ("help", "show help message")
    ("threadcnt,j",
      po::value<int>(&thread_cnt)->default_value(thread_cnt),
      "number of threads")
    ("runcnt",
      po::value<int>(&run_cnt)->default_value(run_cnt),
      "number of runs")
    ("size",
      po::value<int64_t>(&individual_cnt)->default_value(individual_cnt),
      "size of the population")
    ("seed",
      po::value<size_t>(&rand_seed)->default_value(rand_seed),
      "seed for random number generator")
    ("endtime",
      po::value<double>(&end_time)->default_value(end_time),
      "how many days to run")
    ("scenario,s",
      po::value<std::string>(&scenario_file),
      "A NAADSM scenario file")
    ("herd,h",
      po::value<std::string>(&herd_file),
      "A NAADSM XML herd file.")
    ("datafile",
      po::value<std::string>(&data_file)->default_value(data_file),
      "Write to this data file.")
    ("save",
      po::value<bool>(&save_file)->default_value(save_file),
      "Add data to file instead of erasing it with new data.")
    ("loglevel", po::value<std::string>(&log_level)->default_value("info"),
      "Set the logging level to trace, debug, info, warning, error, or fatal.")
    ("translate",
      po::value<std::string>(&translation_file)->default_value(""),
      "write file relating place ids to internal ids")
    ("info", "show provenance of program")
    ;

  for (auto& p : parameters) {
    desc.add_options()(p.name.c_str(),
      po::value<double>(&p.value)->default_value(p.value),
      p.description.c_str());
  }

  po::variables_map vm;
  auto parsed_options=po::parse_command_line(argc, argv, desc);
  po::store(parsed_options, vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  std::map<std::string,std::string> compile_info {
    {"VERSION", VERSION}, {"COMPILETIME", COMPILETIME}
  };
  if (vm.count("info")) {
    for (auto& kv : compile_info) {
      std::cout << kv.second << "\n\n";
    }
    return 0;
  }

  afidd::LogInit(log_level);

  if (test) {
  }

  std::map<ADParam,double*> params;
  for (auto& pm : parameters) {
    params[pm.kind]=&pm.value;
  }

  for (auto& showp : parameters) {
    BOOST_LOG_TRIVIAL(info)<<showp.name<<" "<<showp.value;
  }

  HDFFile file(data_file);
  if (!file.Open(!save_file)) {
    BOOST_LOG_TRIVIAL(error)<<"could not open output file: "<<data_file;
    return -1;
  }
  file.WriteExecutableData(compile_info, parsed_options, {0, 0, 0});
  RandGen rng2(rand_seed+1);

  NAADSMScenario scenario;
  scenario.load(scenario_file, herd_file);
  BOOST_LOG_TRIVIAL(trace) << scenario;

  auto locations=scenario.GetLocations();
  BOOST_LOG_TRIVIAL(debug)<<"locations";
  for ( auto l : locations) {
    double x=std::get<0>(l);
    double y=std::get<1>(l);
    BOOST_LOG_TRIVIAL(debug)<<x<<'\t'<<y;
  }
  file.SaveLocations(locations);

  auto runnable=[=, &scenario](RandGen& rng, size_t single_seed,
      size_t idx)->void {
    std::shared_ptr<TrajectoryObserver> observer=0;
    observer=std::make_shared<TrajectorySave>();

    SIR_run(end_time, parameters, scenario, observer, rng);
    file.SaveTrajectory(parameters, single_seed, idx, observer->Trajectory());
  };

  afidd::smv::Ensemble<decltype(runnable),RandGen> ensemble(runnable,
    thread_cnt, run_cnt, rand_seed);
  ensemble.Run();
  BOOST_LOG_TRIVIAL(debug)<<"Finished running ensemble.";

  file.Close();

  return 0;
}
