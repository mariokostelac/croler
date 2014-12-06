#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "./nucleo_buffer.h"
using std::vector;

int main(int argc, char **argv) {

    vector<char *> string_list;

    FILE *fp = fopen(argv[1], "r");
    if (fp == NULL) {
        perror("Error during opening file.");
        exit(EXIT_FAILURE);
    }
    char *line = NULL;
    size_t len = 0;

    int read_lines = 0;
    while (getline(&line, &len, fp) != -1) {
        char *cpy = new char[strlen(line)];
        snprintf(cpy, strlen(line) - 1, "%s", line);
        string_list.push_back(cpy);
        read_lines++;
    }
    free(line);
    printf("Read %d strings...\n\n", read_lines);
    fclose(fp);

    for (int i = 0, len = string_list.size(); i < len; ++i) {
        NucleoBuffer buff(16);
        for (int j = 0, len = buff.get_size(); j < len; ++j) {
            buff.write(string_list[i][j]);
        }
        printf("%s %04x\n", string_list[i], buff.get_content());
    }

    return 0;
}
