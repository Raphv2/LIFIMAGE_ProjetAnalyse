/***************************************************************************************************
 * @file  analyse.cpp
 * @brief Contains the main program of the project
 **************************************************************************************************/

#include <iostream>
#include "image_io.h"

int main() {
    Image puzzle = read_image("data/puzzle.jpg");
    int width = puzzle.width() * 0.8;

    Color c(0);
    int n = 0;
    for (int i = width ; i < puzzle.width(); ++i) {
        for (int j = 0; j < puzzle.height(); ++j) {
            c = c + srgb(puzzle(i, j));
            ++n;
        }
    }

    c = c/n;

    for (int i = 0 ; i < puzzle.width(); ++i) {
        for (int j = 0; j < puzzle.height(); ++j) {
            puzzle(i, j) = c;
        }
    }

    write_image(puzzle, "data/test.jpg");

    return 0;
}