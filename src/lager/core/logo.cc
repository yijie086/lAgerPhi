// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
// 
// This file is part of lAger.
// 
// lAger is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Shoftware Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// lAger is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with lAger.  If not, see <https://www.gnu.org/licenses/>.
// 

#include "logo.hh"
#include <iostream>

#ifndef LAGER_VERSION
#define LAGER_VERSION "XXX"
#endif

constexpr const size_t LINE_LENGTH = 80;

namespace {
void hline(char c = '_') {
  for (int i = 0; i < LINE_LENGTH; ++i) {
    std::cerr << c;
  }
  std::cerr << "\n\n";
}
void logo1() {
  std::cerr << "______________________________________\n"
               "___  /___    |_  ____/__  ____/__  __ \\\n"
               "__  / __  /| |  / __ __  __/  __  /_/ /\n"
               "_  /___  ___ / /_/ / _  /___  _  _, _/\n"
               "/_____/_/  |_\\____/  /_____/  /_/ |_|\n\n";
}
void logo2() {
  std::cerr
      << "ooooo              .o.         .oooooo.    oooooooooooo ooooooooo.\n"
         "`888'             .888.       d8P'  `Y8b   `888'     `8 `888   "
         "`Y88.\n"
         " 888             .8\"888.     888            888          888   "
         ".d88'\n"
         " 888            .8' `888.    888            888oooo8     888ooo88P'\n"
         " 888           .88ooo8888.   888     ooooo  888    \"     888`88b.\n"
         " 888       o  .8'     `888.  `88.    .88'   888       o  888  `88b.\n"
         "o888ooooood8 o88o     o8888o  `Y8bood8P'   o888ooooood8 o888o  "
         "o888o\n\n";
}

void logo3() {
  std::cerr << "oooo        .o.\n"
               "`888       .888.\n"
               " 888      .8\"888.      .oooooooo  .ooooo.  oooo d8b\n"
               " 888     .8' `888.    888' `88b  d88' `88b `888"
               "8P\n"
               " 888    .88ooo8888.   888   888  888ooo888  888\n"
               " 888   .8'     `888.  `88bod8P'  888    .o  888\n"
               "o888o o88o     o8888o `8oooooo.  `Y8bod8P' d888b\n"
               "                       d\"     YD\n"
               "                       \"Y88888P'\n\n";
}
void info() {
  std::cerr
      << "   Argonne generic l/A-event generator\n"
         "   Version: " LAGER_VERSION "\n"
         "   Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>\n\n"
         "   If you have questions, or want to raise an issue, please go "
         "to:\n"
         "   https://eicweb.phy.anl.gov/monte_carlo/lager\n\n"
         "   If you used lAger to generate data used in a presentation or\n"
         "   an article in a scientific publication, please cite:\n\n"
         "   'S. Joosten, Argonne l/A-event Generator (2021), GitLab\n"
         "    repository, https://eicweb.phy.anl.gov/monte_carlo/lager'\n";
}
}
namespace lager {
void print_logo() {
  hline();
  logo3();
  info();
  hline();
}
}

void logo() { lager::print_logo(); }
