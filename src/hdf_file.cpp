#include <sstream>
#include <cassert>
#include <mutex>
#include <exception>
#include <fstream>
#include "hdf_file.hpp"
#include "hdf5.h"
#include "smv.hpp"


herr_t IterateTrajectories(hid_t group_id, const char* group_name,
  const H5L_info_t* info, void* op_data) {
  int* chosen_idx=static_cast<int*>(op_data);

  std::string sname{group_name};
  if (sname.substr(0,10)==std::string("trajectory")) {
    if (sname.size()>10) {
      std::istringstream convert{sname.substr(10)};
      int val{0};
      try {
        convert>>val;
        *chosen_idx=std::max(*chosen_idx, val+1);
      } catch (std::exception& e) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not convert "<<sname;
      }
    } else {
      *chosen_idx=std::max(*chosen_idx, 1);
    }
  }
  return herr_t{0};
}



class HDFFile::impl {
	hid_t file_id_;
  hid_t trajectory_group_;
  std::string filename_;
  bool open_;
  mutable std::mutex single_writer_;
 public:
  impl(const std::string& filename) : filename_(filename), open_(false),
    file_id_(0), trajectory_group_(0) {}
  ~impl() { this->Close(); }
  bool Open(bool truncate=true) {
    if (truncate) {
      file_id_=H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
          H5P_DEFAULT);
    } else {
      std::ifstream file_exists(filename_.c_str());
      bool exists=file_exists.good();
      file_exists.close();

      if (exists) {
        file_id_=H5Fopen(filename_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        BOOST_LOG_TRIVIAL(debug)<<"File exists: "<<filename_;
      } else {
        BOOST_LOG_TRIVIAL(debug)<<"Creating file: "<<filename_;
        file_id_=H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
      }
    }
    if (file_id_<0) return false;

    open_=true;

    int chosen_idx=0;

    if (!truncate) {
      // Is there already a trajectory in this file?
      herr_t find_status=H5Literate_by_name(file_id_, "/", H5_INDEX_NAME,
        H5_ITER_INC, NULL, IterateTrajectories, &chosen_idx, H5P_DEFAULT);
      BOOST_LOG_TRIVIAL(debug)<<"HDFFile::Open chosen index "<<chosen_idx;
      if (find_status<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not iterate over groups in file.";
      }
    }
    std::stringstream trajname;
    trajname<<"/trajectory";
    if (chosen_idx>0) {
      trajname << chosen_idx;
    }
    BOOST_LOG_TRIVIAL(info)<<"Writing to file directory "<<trajname.str();
    trajectory_group_=H5Gcreate(file_id_, trajname.str().c_str(), H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    return true;
  }


  bool WriteExecutableData(const std::map<std::string,std::string>& compile,
    const boost::program_options::basic_parsed_options<char>& options,
    const std::vector<int64_t>& siri) const {
    std::unique_lock<std::mutex> only_me(single_writer_);

    {
      hsize_t adims=1;
      hid_t dspace_id=H5Screate_simple(1, &adims, NULL);

      for (const auto& kv : compile) {
        hid_t strtype=H5Tcopy(H5T_C_S1);
        herr_t strstatus=H5Tset_size(strtype, kv.second.size());
        if (strstatus<0) {
          BOOST_LOG_TRIVIAL(error)
            <<"Could not create string for executable data.";
            return false;
        }

        hid_t attr0_id=H5Acreate2(trajectory_group_, kv.first.c_str() , strtype,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT);
        herr_t atstatus=H5Awrite(attr0_id, strtype, kv.second.c_str());
        if (atstatus<0) {
          BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<kv.first;
          return false;
        }
        H5Tclose(strtype);
        H5Aclose(attr0_id);
      }
      H5Sclose(dspace_id);
    }

    {
      hsize_t sdims=3;
      hid_t sirspace_id=H5Screate_simple(1, &sdims, NULL);

      hid_t attr1_id=H5Acreate2(trajectory_group_, "Initial Values",
        H5T_STD_I64LE, sirspace_id, H5P_DEFAULT, H5P_DEFAULT);
      herr_t at1status=H5Awrite(attr1_id, H5T_NATIVE_LONG, &siri[0]);
      if (at1status<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute Initial Values";
        return false;
      }
      H5Sclose(sirspace_id);
      H5Aclose(attr1_id);
    }

    // We save the command-line options used to call the program.
    std::stringstream optstring;
    optstring << "<options>";
    for (auto& opt : options.options) {
      optstring << "<option><name>" << opt.string_key << "</name>";
      optstring << "<values>";
      for (auto& v : opt.value) {
        optstring << "<value>" << v << "</value>";
      }
      optstring <<"</values></options>";
    }
    optstring << "</options>";
    
    {
      hsize_t odims=1;
      hid_t ospace_id=H5Screate_simple(1, &odims, NULL);

      hid_t ostrtype=H5Tcopy(H5T_C_S1);
      herr_t ostrstatus=H5Tset_size(ostrtype, optstring.str().size());
      if (ostrstatus<0) {
        BOOST_LOG_TRIVIAL(error)
          <<"Could not create string for executable data.";
          return false;
      }

      hid_t oattr_id=H5Acreate2(trajectory_group_, "Options", ostrtype,
        ospace_id, H5P_DEFAULT, H5P_DEFAULT);
      herr_t ostatus=H5Awrite(oattr_id, ostrtype, optstring.str().c_str());
      if (ostatus<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute Options";
        return false;
      }
      H5Sclose(ospace_id);
      H5Tclose(ostrtype);
      H5Aclose(oattr_id);
    }

    return true;
  }


  bool Close() {
    if (open_) {
      herr_t group_status=H5Gclose(trajectory_group_);
      if (group_status<0) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 group "<<group_status;
        return false;
      }
      herr_t status=H5Fclose(file_id_);
      if (status<0) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 file "<<status;
        return false;
      }
      open_=false;
    }
    return true;
  }


  bool SaveTrajectory(const std::vector<Parameter>& params,
      int seed, int idx, const TrajectoryType& trajectory) const {
    std::unique_lock<std::mutex> only_me(single_writer_);
    assert(open_);

    std::stringstream dset_name;
    dset_name << "dset" << seed << "-" << idx;
    hid_t dset_group=H5Gcreate(trajectory_group_, dset_name.str().c_str(),
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (trajectory.size()==0) {
      BOOST_LOG_TRIVIAL(warning)<<"There was no trajectory to write.";
      H5Gclose(dset_group);
      return false;
    }
    
    std::vector<int64_t> integer(trajectory.size(), 0);

    {
      hsize_t sdims=trajectory.size();
      hid_t sirspace_id=H5Screate_simple(1, &sdims, NULL);
      std::transform(trajectory.begin(), trajectory.end(), integer.begin(),
        [](const TrajectoryEntry& e)->int64_t {
          return e.s;
        });

      hid_t attr1_id=H5Dcreate2(dset_group, "Event",
        H5T_STD_I64LE, sirspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      herr_t at1status=H5Dwrite(attr1_id, H5T_NATIVE_LONG,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, integer.data());
      if (at1status<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute Initial Values";
        return false;
      }
      H5Sclose(sirspace_id);
      H5Dclose(attr1_id);
    }
    {
      hsize_t sdims=trajectory.size();
      hid_t sirspace_id=H5Screate_simple(1, &sdims, NULL);
      std::transform(trajectory.begin(), trajectory.end(), integer.begin(),
        [](const TrajectoryEntry& e)->int64_t {
          return e.i;
        });

      hid_t attr1_id=H5Dcreate2(dset_group, "Whom",
        H5T_STD_I64LE, sirspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      herr_t at1status=H5Dwrite(attr1_id, H5T_NATIVE_LONG,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, integer.data());
      if (at1status<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute Initial Values";
        return false;
      }
      H5Sclose(sirspace_id);
      H5Dclose(attr1_id);
    }
    {
      hsize_t sdims=trajectory.size();
      hid_t sirspace_id=H5Screate_simple(1, &sdims, NULL);
      std::transform(trajectory.begin(), trajectory.end(), integer.begin(),
        [](const TrajectoryEntry& e)->int64_t {
          return e.r;
        });

      hid_t attr1_id=H5Dcreate2(dset_group, "Who",
        H5T_STD_I64LE, sirspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      herr_t at1status=H5Dwrite(attr1_id, H5T_NATIVE_LONG,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, integer.data());
      if (at1status<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute Initial Values";
        return false;
      }
      H5Sclose(sirspace_id);
      H5Dclose(attr1_id);
    }
    {
      hsize_t sdims=trajectory.size();
      hid_t sirspace_id=H5Screate_simple(1, &sdims, NULL);
      std::vector<double> data(trajectory.size(), 0);
      std::transform(trajectory.begin(), trajectory.end(), data.begin(),
        [](const TrajectoryEntry& e)->double {
          return e.t;
        });

      hid_t attr1_id=H5Dcreate2(dset_group, "When",
        H5T_IEEE_F64LE, sirspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      herr_t at1status=H5Dwrite(attr1_id, H5T_NATIVE_DOUBLE,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
      if (at1status<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute Initial Values";
        return false;
      }
      H5Sclose(sirspace_id);
      H5Dclose(attr1_id);
    }

    // Now write dataset attributes.
    hsize_t adims=1;
    hid_t dspace_id=H5Screate_simple(1, &adims, NULL);
    for (auto& p : params) {
      hid_t attr0_id=H5Acreate2(dset_group, p.name.c_str(), H5T_IEEE_F64LE,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      herr_t atstatus=H5Awrite(attr0_id, H5T_NATIVE_DOUBLE, &p.value);
      if (atstatus<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<p.name;
      }
      H5Aclose(attr0_id);
    }
    H5Sclose(dspace_id);

    H5Gclose(dset_group);
    return true;
  }

  bool SaveLocations(const std::vector<std::array<double,2>>& locations) const {
    std::unique_lock<std::mutex> only_me(single_writer_);
    assert(open_);
    hsize_t dims[2];
    dims[0]=locations.size();
    dims[1]=2;
    hid_t dataspace_id=H5Screate_simple(2, dims, NULL);

    hid_t dataset_id=H5Dcreate2(trajectory_group_, "locations",
      H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not create dataset for locations";
      H5Sclose(dataspace_id);
      return false;
    }

    herr_t write_status=H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, locations.data());
    if (write_status<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write data for locations "
        <<write_status;
      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);
      return false;
    }
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return true;
  }

};


HDFFile::HDFFile(const std::string& fname) : pimpl{ new impl{ fname }} {}
HDFFile::~HDFFile() {}
bool HDFFile::Open(bool truncate) { return pimpl->Open(truncate); }
bool HDFFile::Close() { return pimpl->Close(); }
bool HDFFile::SaveTrajectory(const std::vector<Parameter>& params,
  int seed, int idx, const TrajectoryType& traj) const {
  return pimpl->SaveTrajectory(params, seed, idx, traj);
}
bool HDFFile::SaveLocations(const std::vector<std::array<double,2>>& locs) const {
  return pimpl->SaveLocations(locs);
}
bool HDFFile::WriteExecutableData(
    const std::map<std::string,std::string>& compile,
    boost::program_options::basic_parsed_options<char>& cmdline,
    const std::vector<int64_t>& initial_values) const {
  return pimpl->WriteExecutableData(compile, cmdline, initial_values);
}
HDFFile::HDFFile(const HDFFile& o)
: pimpl(o.pimpl) {}
HDFFile& HDFFile::operator=(const HDFFile& o) {
  pimpl=o.pimpl;
  return *this;
}
