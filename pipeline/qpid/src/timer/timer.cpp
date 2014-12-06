/* Copyright 2014 - Mario Kostelac (mario.kostelac@gmail.com) */
#include "./timer.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>

Timer::Timer(const char * n) : name(NULL) {
    int len = strlen(n);
    name = new char[len + 1];
    snprintf(name, len + 1, "%s" , n);
    start_time = clock();
}

Timer::~Timer() {
    delete[] name;
}

Timer* Timer::end(bool print) {
    end_time = clock();

    if (print == true) {
        fprintf(stderr, "* %s finished in %u msec.\n", name,
                (unsigned int) ((end_time - start_time) / (CLOCKS_PER_SEC / 1000.)));
    }

    return this;
}
