#include <vector>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include "./fixed_min_queue.cpp"

int min_elem(std::vector<int>& elems, int size) {
    int min = elems[elems.size() - 1];
    for (int j = 1; j < size; ++j) {
        int curr = elems[elems.size() - j - 1];
        min = curr < min ? curr : min;
    }
    return min;
}

void print_last(std::vector<int>& elems, int size) {
    int sz = elems.size();
    for (int i = size; i >= 1; --i) {
        printf("%d ", elems[sz - i]);
    }
}

void test(int size, int tries, int maxnum) {

    FixedMinQueue<int> pq(size);
    std::vector<int> elems;
    for (int i = 0; i < size; ++i) {
        int num = rand() % (maxnum);
        pq.push(num);
        elems.push_back(num);
    }
    assert(min_elem(elems, size) == pq.min());
    for (int i = 0; i < tries; ++i) {
        int num = rand() % (2 * tries);
        pq.push(num);
        elems.push_back(num);
        assert(min_elem(elems, size) == pq.min());
    }
}

int main(int argc, char** argv) {

    srand(time(NULL));
    for (int i = 0, tests = atoi(argv[1]); i < tests; ++i) {
        test(5, 200, (i % 100) + 1);
        printf("Test %d passed.\n", i);
    }

    return 0;
}
