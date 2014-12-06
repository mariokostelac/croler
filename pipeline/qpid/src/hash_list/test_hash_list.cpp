#include <cassert>
#include <cstdio>
#include "hash_list.cpp"

unsigned int N;
unsigned int KEYS_NUM;

int test_elems_iteration(HashList<int, int>& hl) {
    unsigned int iterated = 0;
    const std::unordered_set<int>& keys = hl.keys();
    for (std::unordered_set<int>::const_iterator key = keys.begin(); key != keys.end(); ++key) {
        HashList<int, int>::List list = hl.get_list(*key);
        HashList<int, int>::List::iterator num = list.begin();
        HashList<int, int>::List::iterator end = list.end();
        for (; num != end; ++num) {
            assert((*num % KEYS_NUM) == *key % KEYS_NUM);
            ++iterated;
        }
    }
    return iterated;
}

void fill_list(HashList<int, int>& list, int size, int keys_num) {
    for (int i = 0; i < size; ++i) {
        list.add(i % keys_num, i);
    }
}

void test_overlapping_keys() {
    HashList<int, int> hl(KEYS_NUM / 2);
    fill_list(hl, N, KEYS_NUM);
    unsigned int iterated = test_elems_iteration(hl);
    assert(iterated == hl.size());
}

int main(int argc, char **argv) {

    if (argc < 3) {
        printf("Callg %s KEYS_NUM ELEMS_NUM\n", argv[0]);
        exit(1);
    }

    HashList<int, int> hl(1000007);
    KEYS_NUM = atoi(argv[1]);
    N = atoi(argv[2]);

    fill_list(hl, N, KEYS_NUM);

    printf("size test: ");
    assert(N == hl.size());
    printf("OK\n");

    printf("key size test: ");
    const std::unordered_set<int>& keys = hl.keys();
    assert(KEYS_NUM == keys.size());
    printf("OK\n");

    printf("elem association test: ");
    unsigned int iterated = test_elems_iteration(hl);
    printf("OK\n");

    printf("elem iterated elems test: ");
    assert(iterated == hl.size());
    printf("OK\n");

    printf("overlapping keys test: ");
    test_overlapping_keys();
    printf("OK\n");

    return 0;
}

