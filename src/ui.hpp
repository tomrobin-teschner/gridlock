#pragma once

#include <string>

#include "ncurses.h"

class UI {
public:
  UI() {
    initscr();
    printw("            _     _ _            _     \n");
    printw("  __ _ _ __(_) __| | | ___   ___| | __ \n");
    printw(" / _` | '__| |/ _` | |/ _ \\ / __| |/ / \n");
    printw("| (_| | |  | | (_| | | (_) | (__|   <  \n");
    printw(" \\__, |_|  |_|\\__,_|_|\\___/ \\___|_|\\_\\ \n");
    printw(" |___/                                 \n");
    refresh();
  }
  ~UI() {
    endwin();
  }
public:
  void draw(int x, int y, std::string text) {
    mvprintw(y, x, text.c_str());
    refresh();
  }
  void block() {
    getch();
  }
};