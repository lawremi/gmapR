/* Functions for reading data from byte streams. Inlined because they
   are tiny and need to be fast. */

static inline int
read_int (unsigned char **bytes) {
  int x = ((int *)*bytes)[0];
  *bytes += sizeof(int);
  return x;
}

static inline char
read_char (unsigned char **bytes) {
  char x = (*bytes)[0];
  (*bytes)++;
  return x;
}

static inline const char *
read_string (unsigned char **bytes) {
  const char *string = *bytes;
  (*bytes) += strlen(string) + 1;
  return string;
}
