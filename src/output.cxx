// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "output.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options/variables_map.hpp>

#include "event.h"

namespace trento {

namespace {

void write_stream(std::ostream& os,
    int num, double impact_param, const Event& event) {
  using std::fixed;
  using std::setprecision;
  using std::setw;
  using std::scientific;

  os << setprecision(10)       << num
     << setw(15) << fixed      << impact_param
     << setw(5)                << event.npart()
     << setw(18) << scientific << event.multiplicity()
     << fixed;

  for (const auto& ecc : event.eccentricity())
    os << setw(14)             << ecc.second;

  os << '\n';
}

void write_text_file(const fs::path& path,
    int num, double impact_param, const Event& event) {
  fs::ofstream ofs{path};

  // Write a commented header of event properties as key = value pairs.
  ofs << std::setprecision(10)
      << "# event "   << num                  << '\n'
      << "# b     = " << impact_param         << '\n'
      << "# npart = " << event.npart()        << '\n'
      << "# mult  = " << event.multiplicity() << '\n';

  for (const auto& ecc : event.eccentricity())
    ofs << "# e" << ecc.first << "    = " << ecc.second << '\n';

  // Write IC profile as a block grid.  Use C++ default float format (not
  // fixed-width) so that trailing zeros are omitted.  This significantly
  // increases output speed and saves disk space since many grid elements are
  // zero.
  for (const auto& row : event.reduced_thickness_grid()) {
    auto&& iter = row.begin();
    // Write all row elements except the last with a space delimiter afterwards.
    do {
      ofs << *iter << ' ';
    } while (++iter != --row.end());

    // Write the last element and a linebreak.
    ofs << *iter << '\n';
  }
}

}  // unnamed namespace

Output::Output(const VarMap& var_map) {
  auto nevents = var_map["number-events"].as<int>();
  auto width = static_cast<int>(std::ceil(std::log10(nevents)));

  // Write to stdout unless the quiet option was specified.
  if (!var_map["quiet"].as<bool>()) {
    writers_.emplace_back(
      [width](int num, double impact_param, const Event& event) {
        std::cout << std::setw(width);
        write_stream(std::cout, num, impact_param, event);
      }
    );
  }

  if (var_map.count("output")) {
    auto output_path = var_map["output"].as<fs::path>();
    if (output_path.has_extension() && (
          output_path.extension() == ".hdf5" ||
          output_path.extension() == ".hd5"  ||
          output_path.extension() == ".h5"
          )) {
      throw std::runtime_error{"HDF5 not implemented yet"};  // TODO
    } else {
      if (fs::exists(output_path)) {
        if (!fs::is_empty(output_path)) {
          throw std::runtime_error{"output directory '" + output_path.string() +
                                   "' must be empty"};
        }
      } else {
        fs::create_directories(output_path);
      }
      writers_.emplace_back(
        [output_path, width](int num, double impact_param, const Event& event) {
          std::ostringstream padded_fname{};
          padded_fname << std::setw(width) << std::setfill('0')
                       << num << ".dat";
          write_text_file(output_path / padded_fname.str(),
              num, impact_param, event);
        }
      );
    }
  }
}

}  // namespace trento
