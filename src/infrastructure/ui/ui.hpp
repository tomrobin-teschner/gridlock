#pragma once

#include <string>

#include "ncurses.h"

#include "src/infrastructure/numToStringConverter/numToStringConverter.hpp"

class UI {
public:
  UI(int maxTimeSteps, int maxPicardIterations);
  ~UI();

public:
  void createSkeleton();
  void updateTimestatistics(int time, double CFL, double dt, double totalTime, int ehh, int emm, int ess,
    int rhh, int rmm, int rss);
  void updatePicardIteration(int picardIterations, double uRes, double vRes, double pRes, int uIter, int vIter,
    int pIter);
  void updateIteration(double uRes, double vRes, double pRes, double divU);
  void draw(int x, int y, std::string text);
  void block();
private:
  int _maxTimeSteps;
  int _maxPicardIterations;
};