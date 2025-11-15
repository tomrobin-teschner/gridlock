#include "src/infrastructure/ui/ui.hpp"

UI::UI(Parameters params) : _maxTimeSteps(params.solver<int>("time", "timeSteps")),
  _maxPicardIterations(params.solver<int>("linearisation", "maxPicardIterations")) {

  initscr();
  printw("            _     _ _            _     \n");
  printw("  __ _ _ __(_) __| | | ___   ___| | __ \n");
  printw(" / _` | '__| |/ _` | |/ _ \\ / __| |/ / \n");
  printw("| (_| | |  | | (_| | | (_) | (__|   <  \n");
  printw(" \\__, |_|  |_|\\__,_|_|\\___/ \\___|_|\\_\\ \n");
  printw(" |___/                                 \n");

  this->draw(2, 42, "Gridlock 0.11.0");
  this->draw(3, 42, "A 2D structured CFD solver for incompressible flows");
  this->draw(4, 42, "Tom-Robin Teschner");

  refresh();
}

UI::~UI() {
  endwin();
}

void UI::createSkeleton() {
  // time information
  mvprintw(8, 0, "Time step");
  mvprintw(8, 15, "CFL");
  mvprintw(8, 30, "dt (s)");
  mvprintw(8, 45, "Total time (s)");
  mvprintw(8, 60, "Elapsed");
  mvprintw(8, 75, "Remaining");

  // picard iteration information
  mvprintw(11, 0, "Picard loop");
  mvprintw(11, 15, "u-residual");
  mvprintw(11, 30, "v-residual");
  mvprintw(11, 45, "p-residual");
  mvprintw(11, 60, "u-iterations");
  mvprintw(11, 75, "v-iterations");
  mvprintw(11, 90, "p-iterations");

  // end of timestep residuals
  mvprintw(14, 0, "u-residual");
  mvprintw(14, 15, "v-residual");
  mvprintw(14, 30, "p-residual");
  mvprintw(14, 45, "div(u)");

  refresh();
}

void UI::updateTimestatistics(int time, double CFL, double dt, double totalTime, int ehh, int emm, int ess,
  int rhh, int rmm, int rss) {
  this->draw(9, 0, NumToStringConverter::integer(time) + " / " + NumToStringConverter::integer(_maxTimeSteps));
  this->draw(9, 15, NumToStringConverter::fixedFloat(CFL, 1));
  this->draw(9, 30, NumToStringConverter::scientificFloat(dt, 2));
  this->draw(9, 45, NumToStringConverter::scientificFloat(totalTime, 2));
  this->draw(9, 60, NumToStringConverter::time(ehh, emm, ess));
  this->draw(9, 75, NumToStringConverter::time(rhh, rmm, rss));
}

void UI::updatePicardIteration(int picardIterations, double uRes, double vRes, double pRes, int uIter, int vIter,
  int pIter) {
  this->draw(12, 0, NumToStringConverter::integer(picardIterations) + " / " +
    NumToStringConverter::integer(_maxPicardIterations));
  this->draw(12, 15, NumToStringConverter::scientificFloat(uRes, 2));
  this->draw(12, 30, NumToStringConverter::scientificFloat(vRes, 2));
  this->draw(12, 45, NumToStringConverter::scientificFloat(pRes, 2));
  this->draw(12, 60, NumToStringConverter::integer(uIter));
  this->draw(12, 75, NumToStringConverter::integer(vIter));
  this->draw(12, 90, NumToStringConverter::integer(pIter));
}

void UI::updateIteration(double uRes, double vRes, double pRes, double divU) {
  this->draw(15, 0, NumToStringConverter::scientificFloat(uRes, 2));
  this->draw(15, 15, NumToStringConverter::scientificFloat(vRes, 2));
  this->draw(15, 30, NumToStringConverter::scientificFloat(pRes, 2));
  this->draw(15, 45, NumToStringConverter::scientificFloat(divU, 2));
}

void UI::draw(int x, int y, std::string text) {
  std::string blanks(15, ' ');
  mvprintw(x, y, blanks.c_str());
  mvprintw(x, y, text.c_str());
  refresh();
}

void UI::block() {
  getch();
}