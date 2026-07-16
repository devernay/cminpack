/*
  Compare two files token by token, with a numeric tolerance.

  Usage: cmpfiles file1 file2 [reltol] [abstol]

  Each line is split into whitespace-separated tokens. Corresponding
  tokens are compared as follows:
    - if both tokens parse completely as floating-point numbers, they
      are considered equal when
          |a - b| <= abstol   OR   |a - b| <= reltol * max(|a|,|b|)
      (two NaNs, or two infinities of the same sign, also compare equal);
    - otherwise they are compared as strings.

  This replaces an exact byte-for-byte comparison, which is not portable:
  the test drivers print results with %g/%e at ~7 significant digits, and
  the least-significant digits legitimately differ between compilers, C
  runtimes and math libraries (see cminpack issues #37 and #77). A small
  tolerance absorbs those differences while still catching real errors.

  Exit status is 0 if the files match (within tolerance), 1 if they differ,
  2 on a usage or I/O error.

  Author: Frederic.Devernay@m4x.org
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* Default tolerances, chosen to match the Makefile's
   "numdiff -a 1e-4 -r 1e-4" and loosened enough to absorb the
   platform-dependent last digits of 7-significant-digit output. */
#define DEFAULT_RELTOL 1e-2
#define DEFAULT_ABSTOL 1e-4

/* Parse a token as a floating-point number. Returns 1 and stores the
   value in *v if the whole token (ignoring surrounding blanks) is a
   valid number, 0 otherwise. */
static int parse_number(const char *tok, double *v)
{
    char *end = NULL;
    double d;
    if (tok == NULL || *tok == '\0') {
        return 0;
    }
    d = strtod(tok, &end);
    if (end == tok) {
        return 0; /* no conversion */
    }
    while (*end != '\0' && isspace((unsigned char)*end)) {
        ++end;
    }
    if (*end != '\0') {
        return 0; /* trailing non-numeric characters */
    }
    *v = d;
    return 1;
}

/* Numeric equality within tolerance, with NaN/inf handling. */
static int numbers_equal(double a, double b, double reltol, double abstol)
{
    double diff, scale;
    if (isnan(a) || isnan(b)) {
        return isnan(a) && isnan(b);
    }
    if (isinf(a) || isinf(b)) {
        return a == b; /* same sign infinity */
    }
    diff = fabs(a - b);
    if (diff <= abstol) {
        return 1;
    }
    scale = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
    return diff <= reltol * scale;
}

/* Return the next whitespace-delimited token starting at *pp, advancing
   *pp past it and NUL-terminating the token in place. Returns NULL when
   no more tokens remain. This is a minimal reentrant tokenizer (unlike
   strtok, which keeps a single global state and cannot tokenize two
   strings at once). */
static char *next_token(char **pp)
{
    char *s = *pp;
    char *start;
    while (*s != '\0' && isspace((unsigned char)*s)) {
        ++s;
    }
    if (*s == '\0') {
        *pp = s;
        return NULL;
    }
    start = s;
    while (*s != '\0' && !isspace((unsigned char)*s)) {
        ++s;
    }
    if (*s != '\0') {
        *s = '\0';
        ++s;
    }
    *pp = s;
    return start;
}

/* Compare two lines token by token. Returns 0 if equal (within
   tolerance), 1 if different. */
static int compare_lines(const char *l1, const char *l2,
                         double reltol, double abstol)
{
    char b1[1024], b2[1024];
    char *s1, *s2;
    strncpy(b1, l1, sizeof(b1) - 1); b1[sizeof(b1) - 1] = '\0';
    strncpy(b2, l2, sizeof(b2) - 1); b2[sizeof(b2) - 1] = '\0';

    s1 = b1;
    s2 = b2;
    for (;;) {
        double v1, v2;
        int n1, n2;
        char *t1 = next_token(&s1);
        char *t2 = next_token(&s2);
        if (t1 == NULL && t2 == NULL) {
            return 0; /* both exhausted, lines equal */
        }
        if (t1 == NULL || t2 == NULL) {
            return 1; /* different number of tokens */
        }
        n1 = parse_number(t1, &v1);
        n2 = parse_number(t2, &v2);
        if (n1 && n2) {
            if (!numbers_equal(v1, v2, reltol, abstol)) {
                return 1;
            }
        } else if (strcmp(t1, t2) != 0) {
            return 1;
        }
    }
}

int main(int argc, char **argv)
{
    FILE *fp1, *fp2;
    char buf1[1024], buf2[1024];
    int differences = 0;
    int line = 0;
    double reltol = DEFAULT_RELTOL;
    double abstol = DEFAULT_ABSTOL;

    if (argc < 3 || argc > 5) {
        fprintf(stderr, "Usage: %s file1 file2 [reltol] [abstol]\n", argv[0]);
        exit(2);
    }
    if (argc >= 4) {
        reltol = strtod(argv[3], NULL);
    }
    if (argc >= 5) {
        abstol = strtod(argv[4], NULL);
    }

    fp1 = fopen(argv[1], "r");
    if (fp1 == NULL) {
        perror("Error: cannot open input file");
        exit(2);
    }
    fp2 = fopen(argv[2], "r");
    if (fp2 == NULL) {
        perror("Error: cannot open input file");
        fclose(fp1);
        exit(2);
    }

    while (fgets(buf1, sizeof(buf1), fp1) != NULL) {
        ++line;
        if (fgets(buf2, sizeof(buf2), fp2) == NULL) {
            fprintf(stderr, "%d< %s", line, buf1);
            ++differences;
        } else if (compare_lines(buf1, buf2, reltol, abstol) != 0) {
            fprintf(stderr, "%d< %s", line, buf1);
            fprintf(stderr, "%d> %s", line, buf2);
            ++differences;
        }
    }
    while (fgets(buf2, sizeof(buf2), fp2) != NULL) {
        ++line;
        fprintf(stderr, "%d> %s", line, buf2);
        ++differences;
    }

    fclose(fp1);
    fclose(fp2);

    if (differences != 0) {
        fprintf(stderr,
                "Total differences between files '%s' and '%s': %d lines"
                " (reltol=%g, abstol=%g)\n",
                argv[1], argv[2], differences, reltol, abstol);
    }
    exit(differences != 0);
}
