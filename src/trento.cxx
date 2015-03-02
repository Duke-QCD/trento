// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/any.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "config.h"
#include "fwd_decl.h"

// CMake sets this definition.
// Fall back to a sane default.
#ifndef TRENTO_VERSION_STRING
#define TRENTO_VERSION_STRING "dev"
#endif

namespace trento {

namespace {

void print_short_usage(std::ostream& os = std::cout) {
  os << "usage: trento [options] projectile projectile [number-events = 1]\n";
}

void print_long_usage() {
  print_short_usage();
  std::cout <<
    "\n"
    "Run the trento initial condition model.\n"
    "\n"
    "Projectiles must be standard symbols, e.g. Pb for lead and p for proton.\n"
    "This would run one thousand proton-lead events with default settings:\n"
    "\n"
    "  trento p Pb 1000\n"
    "\n"
    "The default output mode is to print these event properties to stdout:\n"
    "\n"
    "  event_number total_entropy impact_param Npart e_2 e_3 e_4\n"
    "\n"
    "where the e_n are eccentricity harmonics. Give the -q/--quiet option to disable\n"
    "this output.\n"
    "\n"
    "Event profiles are not output by default; the -o/--output option must be\n"
    "explicitly given. If an HDF5 filename is given (.h5, .hd5, .hdf5), all profiles\n"
    "will be written to file. Otherwise, the option is interpreted as a directory and\n"
    "profiles will be written to numbered text files in the directory.\n"
    "\n"
    "All options (including projectiles and number of events) can be saved in\n"
    "configuration files and used via the -c/--config-file option. If an option is\n"
    "specified in both a config file and on the command line, the command line\n"
    "overrides the config file.\n"
    "\n"
    "See the online documentation for complete usage information.\n";
}

void print_version() {
  std::cout << "trento " << TRENTO_VERSION_STRING << '\n';
}

void print_bibtex() {
  std::cout <<
    "@article{Moreland:2014oya,\n"
    "      author         = \"Moreland, J. Scott and Bernhard, Jonah E. and Bass,\n"
    "                        Steffen A.\",\n"
    "      title          = \"{An effective model for entropy deposition in high-energy\n"
    "                        pp, pA, and AA collisions}\",\n"
    "      year           = \"2014\",\n"
    "      eprint         = \"1412.4708\",\n"
    "      archivePrefix  = \"arXiv\",\n"
    "      primaryClass   = \"nucl-th\",\n"
    "      SLACcitation   = \"%%CITATION = ARXIV:1412.4708;%%\",\n"
    "}\n";
}

void print_default_config() {
  std::cout <<
    "# default trento configuration\n"
    "projectile = \n"
    "projectile = \n"
    "number-events = 1\n"
    "\n"
    "quiet = false\n"
    "# output = \n"
    "\n";
}

}  // unnamed namespace

}  // namespace trento

int main(int argc, char* argv[]) {
  using namespace trento;

  // Parse options with boost::program_options.
  // There are quite a few options, so let's separate them into logical groups.
  using opt_desc = po::options_description;

  using vec_str = std::vector<std::string>;
  opt_desc main_opts{};
  main_opts.add_options()
    ("projectile", po::value<vec_str>()->required(),
     "projectile symbols")
    ("number-events", po::value<std::size_t>()->default_value(1),
     "number of events");

  // Make all main arguments positional.
  po::positional_options_description positional_opts{};
  positional_opts
    .add("projectile", 2)
    .add("number-events", 1);

  using vec_path = std::vector<fs::path>;
  opt_desc general_opts{"general options"};
  general_opts.add_options()
    ("help,h", "show this help message and exit")
    ("version", "print version information and exit")
    ("bibtex", "print bibtex entry and exit")
    ("default-config", "print a config file with default settings and exit")
    ("config-file,c", po::value<vec_path>()->value_name("FILE"),
     "configuration file, can be passed multiple times");

  opt_desc output_opts{"output options"};
  output_opts.add_options()
    ("quiet,q", po::bool_switch(),
     "do not print event properties to stdout")
    ("output,o", po::value<fs::path>()->value_name("PATH"),
     "HDF5 file or directory for text files");

  opt_desc phys_opts{"physical options"};
  phys_opts.add_options()
    ("reduced-thickness,p",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "reduced thickness parameter")
    ("fluctuation,k",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
     "gamma fluctuation shape parameter")
    ("beam-energy,s",
     po::value<double>()->value_name("FLOAT")->default_value(2760., "2760"),
     "beam energy sqrt(s) [GeV]")
    ("nucleon-width,w",
     po::value<double>()->value_name("FLOAT")->default_value(.6, "0.6"),
     "Gaussian nucleon width [fm]")
    ("cross-section,x",
     po::value<double>()->value_name("FLOAT")->default_value(6.4, "6.4"),
     "inelastic nucleon-nucleon cross section sigma_NN [fm^2]")
    ("normalization,n",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
     "normalization factor")
    ("b-min,b",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum impact parameter [fm]")
    ("b-max,B",
     po::value<double>()->value_name("FLOAT")->default_value(0., "auto"),
     "maximum impact parameter [fm]");

  opt_desc grid_opts{"grid options"};
  grid_opts.add_options()
    ("grid-size",
     po::value<double>()->value_name("FLOAT")->default_value(10., "10.0"),
     "grid size [fm], grid will extend from -size to +size")
    ("grid-step",
     po::value<double>()->value_name("FLOAT")->default_value(.1, "0.1"),
     "grid step size [fm]");

  // Make a meta-group containing all the option groups except the main
  // positional options (don't want the auto-generated usage info for those).
  opt_desc usage_opts{};
  usage_opts
    .add(general_opts)
    .add(output_opts)
    .add(phys_opts)
    .add(grid_opts);

  // Initialize this before the try...catch block so it stays in scope
  // afterwards.
  po::variables_map var_map;

  try {
    // Now make a meta-group containing _all_ options.
    opt_desc all_opts{};
    all_opts
      .add(usage_opts)
      .add(main_opts);

    // Parse command line options.
    po::store(po::command_line_parser(argc, argv)
        .options(all_opts).positional(positional_opts).run(), var_map);

    //
    if (var_map.count("help")) {
      print_long_usage();
      std::cout << usage_opts;
      return 0;
    }
    if (var_map.count("version")) {
      print_version();
      return 0;
    }
    if (var_map.count("bibtex")) {
      print_bibtex();
      return 0;
    }
    if (var_map.count("default-config")) {
      print_default_config();
      return 0;
    }

    if (var_map.count("config-file")) {
      opt_desc config_file_opts{};
      config_file_opts
        .add(main_opts)
        .add(output_opts)
        .add(phys_opts)
        .add(grid_opts);

      for (const auto &path : var_map["config-file"].as<vec_path>()) {
        if (! fs::exists(path)) {
          std::cerr << "configuration file " << path << " not found\n";
          return 1;
        }
        fs::ifstream ifs{path};
        po::store(po::parse_config_file(ifs, config_file_opts), var_map);
      }
    }

    //
    po::notify(var_map);

    if (var_map["projectile"].as<vec_str>().size() != 2)
      throw po::required_option{"projectile"};
  }
  catch (const po::required_option& e) {
    //
    print_short_usage(std::cerr);
    std::cerr <<
      "\n"
      "must specify exactly two projectiles\n"
      "\n"
      "run 'trento --help' for more information\n";
    return 1;
  }
  catch (const std::exception& e) {
    //
    print_short_usage(std::cerr);
    std::cerr
      << '\n' << e.what()
      << "\n\nrun 'trento --help' for more information\n";
    return 1;
  }

  for (const auto& key : {
      "fluctuation",
      "beam-energy",
      "nucleon-width",
      "cross-section",
      "normalization",
      "grid-size",
      "grid-step",
      }) {
    if (var_map[key].as<double>() <= 0.) {
      std::cerr << key << " must be positive\n";
      return 1;
    }
  }

  for (const auto& key : {"b-min", "b-max"}) {
    if (var_map[key].as<double>() < 0.) {
      std::cerr << key << " must be non-negative\n";
      return 1;
    }
  }

  auto projectiles = var_map["projectile"].as<vec_str>();
  std::cout
    << "projectiles   = "
    << projectiles[0] << ' '
    << projectiles[1] << '\n'
    << "number events = " << var_map["number-events"].as<std::size_t>() << '\n'
    << '\n'
    << "quiet = " << std::boolalpha << var_map["quiet"].as<bool>() << '\n';

  if (var_map.count("output"))
    std::cout << "output = " << var_map["output"].as<fs::path>() << '\n';

  std::cout << '\n';

  for (const auto& opt : phys_opts.options()) {
    auto key = opt->key("*");
    std::cout << key << " = " << var_map[key].as<double>() << '\n';
  }

  std::cout << '\n';

  for (const auto& opt : grid_opts.options()) {
    auto key = opt->key("*");
    std::cout << key << " = " << var_map[key].as<double>() << '\n';
  }

  set_static_vars(var_map);

  return 0;
}
