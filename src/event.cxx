// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "event.h"

#include <iomanip>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options/variables_map.hpp>

namespace trento {

Event::Event(int grid_steps)
    : TA(boost::extents[grid_steps][grid_steps]),
      TB(boost::extents[grid_steps][grid_steps]),
      TR(boost::extents[grid_steps][grid_steps])
{}

namespace {

std::ostream& operator<<(std::ostream& os, const Event& event) {
  using std::fixed;
  using std::setprecision;
  using std::setw;
  using std::scientific;

  os << setprecision(10)       << event.num
     << setw(15) << fixed      << event.impact_param
     << setw(5)                << event.npart
     << setw(18) << scientific << event.multiplicity
     << fixed;

  for (const auto& ecc : event.eccentricity)
    os << setw(14)             << ecc.second;

  return os << '\n';
}

void write_text_file(const fs::path& path, const Event& event) {
  fs::ofstream ofs{path};

  // Write a commented header of event properties as key = value pairs.
  ofs << std::setprecision(10)
      << "# event "   << event.num          << '\n'
      << "# b     = " << event.impact_param << '\n'
      << "# npart = " << event.npart        << '\n'
      << "# mult  = " << event.multiplicity << '\n';

  for (const auto& ecc : event.eccentricity)
    ofs << "# e" << ecc.first << "    = " << ecc.second << '\n';

  // Write IC profile as a block grid.  Use C++ default float format (not
  // fixed-width) so that trailing zeros are omitted.  This significantly
  // increases output speed and saves disk space since many grid elements are
  // zero.
  for (const auto& row : event.TR) {
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

OutputFunctionVector create_output_functions(const VarMap& var_map) {
  OutputFunctionVector output_functions;

  auto num_events = var_map["number-events"].as<int>();
  auto event_num_width = static_cast<int>(std::ceil(std::log10(num_events)));

  // Write to stdout unless the quiet option was specified.
  if (!var_map["quiet"].as<bool>()) {
    output_functions.emplace_back(
      [event_num_width](const Event& event) {
        std::cout << std::setw(event_num_width) << event;
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
      output_functions.emplace_back(
        [output_path, event_num_width](const Event& event) {
          std::ostringstream padded_fname{};
          padded_fname << std::setw(event_num_width) << std::setfill('0')
                       << event.num << ".dat";
          write_text_file(output_path / padded_fname.str(), event);
        }
      );
    }
  }

  return output_functions;
}

}  // namespace trento
