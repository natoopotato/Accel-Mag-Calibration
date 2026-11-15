#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathlib.h"
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <stdarg.h>

#define MIN_POINTS 10
#define USAGE "Usage: %s -m <magnitude> -f <filename> [-o] [-j] -h"

/* Outlier filter settings (-o). Both use robust statistics (median and median
 * absolute deviation), because a mean and standard deviation are themselves
 * pulled around by the bad points. 1.4826 scales a MAD to a standard deviation. */
#define MAD_TO_SIGMA 1.4826
#define PREFILTER_SIGMAS 5.0     /* coarse pass on raw data, catches huge values */
#define RESIDUAL_SIGMAS 3.5      /* fine pass on the fit residuals */
#define RESIDUAL_FLOOR 0.005     /* never reject within 0.5% of the target magnitude */
#define MAX_FILTER_PASSES 5

/* Read one line of any length into *buf (grown with realloc), without the
 * trailing '\n'. Returns 1 if a line was read, 0 on EOF with nothing read,
 * -1 on allocation failure. */
static int read_line(FILE *fp, char **buf, size_t *cap, size_t *len) {
    int c;
    *len = 0;
    while ((c = getc(fp)) != EOF) {
        if (*len + 1 >= *cap) {
            size_t ncap = *cap ? *cap * 2 : 256;
            char *nb = (char *) realloc(*buf, ncap);
            if (nb == NULL)
                return -1;
            *buf = nb;
            *cap = ncap;
        }
        if (c == '\n') {
            (*buf)[*len] = '\0';
            return 1;
        }
        (*buf)[(*len)++] = (char) c;
    }
    if (*len == 0)
        return 0;
    (*buf)[*len] = '\0';
    return 1;
}

static const char *skip_blanks(const char *p) {
    while (*p == ' ' || *p == '\t')
        p++;
    return p;
}

/* Parse "x y z" separated by spaces/tabs with an optional single comma.
 * Returns 1 on success, 0 on any malformed or non-finite input. */
static int parse_point(const char *s, double *v) {
    const char *p = skip_blanks(s);
    int k;
    for (k = 0; k < 3; k++) {
        char *end;
        if (k > 0) {
            const char *q = skip_blanks(p);
            if (*q == ',')
                q = skip_blanks(q + 1);
            if (q == p)
                return 0; /* no separator */
            p = q;
        }
        if (*p == '\0')
            return 0;
        v[k] = strtod(p, &end);
        if (end == p || !isfinite(v[k]))
            return 0;
        p = end;
    }
    while (*p != '\0') {
        if (!isspace((unsigned char) *p))
            return 0;
        p++;
    }
    return 1;
}

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *) a, y = *(const double *) b;
    return (x > y) - (x < y);
}

/* Median of a[0..n-1]. Sorts a in place. */
static double median_sorted(double *a, int n) {
    qsort(a, (size_t) n, sizeof(double), cmp_double);
    if (n % 2)
        return a[n / 2];
    return 0.5 * (a[n / 2 - 1] + a[n / 2]);
}

/* Keep only the points with keep[i] set, compacting pts in place.
 * Returns the new point count. */
static int compact_points(double *pts, int n, const char *keep) {
    int i, m = 0;
    for (i = 0; i < n; i++) {
        if (keep[i]) {
            if (m != i) {
                pts[3 * m] = pts[3 * i];
                pts[3 * m + 1] = pts[3 * i + 1];
                pts[3 * m + 2] = pts[3 * i + 2];
            }
            m++;
        }
    }
    return m;
}

/* Coarse pass on the raw data, before any fit. A huge value wrecks the first
 * ellipsoid fit (the fit squares every coordinate), so the residual filter below
 * cannot be trusted until gross errors are gone. Estimates the centre as the
 * per-axis median, then drops points whose distance from that centre is far
 * above the median distance. Only the high side is tested: a point that is too
 * close to the centre is not a gross error. Returns 0 on success (n updated),
 * -1 on allocation failure. */
static int prefilter_outliers(double *pts, int *n) {
    int i, k, kept = 0;
    double *tmp = (double *) malloc((size_t) *n * sizeof(double));
    double *rad = (double *) malloc((size_t) *n * sizeof(double));
    char *keep = (char *) malloc((size_t) *n);
    double c[3], med, mad, limit;

    if (tmp == NULL || rad == NULL || keep == NULL) {
        free(tmp);
        free(rad);
        free(keep);
        return -1;
    }
    for (k = 0; k < 3; k++) {
        for (i = 0; i < *n; i++)
            tmp[i] = pts[3 * i + k];
        c[k] = median_sorted(tmp, *n);
    }
    for (i = 0; i < *n; i++) {
        double dx = pts[3 * i] - c[0], dy = pts[3 * i + 1] - c[1], dz = pts[3 * i + 2] - c[2];
        rad[i] = sqrt(dx * dx + dy * dy + dz * dz);
        tmp[i] = rad[i];
    }
    med = median_sorted(tmp, *n);
    for (i = 0; i < *n; i++)
        tmp[i] = fabs(rad[i] - med);
    mad = median_sorted(tmp, *n);
    limit = med + PREFILTER_SIGMAS * MAD_TO_SIGMA * mad;

    if (mad > 0.0) {
        for (i = 0; i < *n; i++) {
            keep[i] = (char) (rad[i] <= limit);
            kept += keep[i];
        }
        if (kept >= MIN_POINTS && kept < *n)
            *n = compact_points(pts, *n, keep);
        else if (kept < MIN_POINTS)
            fprintf(stderr, "Outlier filter skipped: it would leave fewer than %d points.\n", MIN_POINTS);
    }
    free(tmp);
    free(rad);
    free(keep);
    return 0;
}

/* Fine pass, after a fit. Calibrates every point with (B, A_1) and measures how far
 * its length is from the target magnitude. Points whose error is far from the
 * median error are dropped. Returns the number of points dropped (n updated),
 * or -1 on allocation failure. */
static int residual_filter(double *pts, int *n, const double *B, const double *A_1, double magnitude) {
    int i, kept = 0, dropped;
    double *res = (double *) malloc((size_t) *n * sizeof(double));
    double *tmp = (double *) malloc((size_t) *n * sizeof(double));
    char *keep = (char *) malloc((size_t) *n);
    double med, mad, limit;

    if (res == NULL || tmp == NULL || keep == NULL) {
        free(res);
        free(tmp);
        free(keep);
        return -1;
    }
    for (i = 0; i < *n; i++) {
        double d0 = pts[3 * i] - B[0], d1 = pts[3 * i + 1] - B[1], d2 = pts[3 * i + 2] - B[2];
        double c0 = A_1[0] * d0 + A_1[1] * d1 + A_1[2] * d2;
        double c1 = A_1[3] * d0 + A_1[4] * d1 + A_1[5] * d2;
        double c2 = A_1[6] * d0 + A_1[7] * d1 + A_1[8] * d2;
        res[i] = sqrt(c0 * c0 + c1 * c1 + c2 * c2) - magnitude;
        tmp[i] = res[i];
    }
    med = median_sorted(tmp, *n);
    for (i = 0; i < *n; i++)
        tmp[i] = fabs(res[i] - med);
    mad = median_sorted(tmp, *n);
    limit = RESIDUAL_SIGMAS * MAD_TO_SIGMA * mad;
    if (limit < RESIDUAL_FLOOR * magnitude)
        limit = RESIDUAL_FLOOR * magnitude;

    for (i = 0; i < *n; i++) {
        keep[i] = (char) (fabs(res[i] - med) <= limit);
        kept += keep[i];
    }
    dropped = *n - kept;
    if (dropped > 0 && kept < MIN_POINTS) {
        fprintf(stderr, "Outlier filter stopped: it would leave fewer than %d points.\n", MIN_POINTS);
        dropped = 0;
    } else if (dropped > 0) {
        *n = compact_points(pts, *n, keep);
    }
    free(res);
    free(tmp);
    free(keep);
    return dropped;
}

/* Print s as a JSON string literal (quotes and escapes included). */
static void print_json_string(const char *s) {
    putchar('"');
    for (; *s != '\0'; s++) {
        unsigned char c = (unsigned char) *s;
        if (c == '"' || c == '\\')
            printf("\\%c", c);
        else if (c < 0x20)
            printf("\\u%04x", c);
        else
            putchar(c);
    }
    putchar('"');
}

static int fit_fail(int quiet, const char *fmt, ...) {
    if (!quiet) {
        va_list ap;
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        va_end(ap);
    }
    return 1;
}

/* Fit an ellipsoid to the n points in pts (x y z each) and return the hard-iron
 * offset B[3] and the soft-iron matrix A_1[9]. Returns 0 on success. On failure
 * returns 1 and, unless quiet, prints the reason to stderr. */
static int fit_ellipsoid(const double *pts, int n, double magnitude, double *B_out, double *A_1_out, int quiet) {
    double *D;
    double S[100], S11[36], S12[24], S12t[24], S22[16], S22_1[16], S22a[24], S22b[36];
    double SS[36], E[36], SSS[36], U[3];
    /* Pre-inverted constraint matrix */
    double C[36] = {
        0.0, 0.5, 0.5, 0.0, 0.0, 0.0,
        0.5, 0.0, 0.5, 0.0, 0.0, 0.0,
        0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, -0.25, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, -0.25, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, -0.25
    };
    double eigen_real[6], eigen_imag[6], v1[6], v2[4], v[10], Q[9], Q_1[9], B[3], QB[3];
    double J, hmb, SSSS[9];
    int i, index;
    double maxval, norm, btqb, eigen_real3[3], eigen_imag3[3], Dz[9], vdz[9], SQ[9], A_1[9];
    double norm1, norm2, norm3;
    double x, y, z;

    if ((size_t) n > SIZE_MAX / (10 * sizeof(double)))
        return fit_fail(quiet, "too many data points\n");
    D = (double *) malloc(10 * (size_t) n * sizeof(double));
    if (D == NULL)
        return fit_fail(quiet, "Out of memory\n");
    for (i = 0; i < n; i++) {
        x = pts[3 * (size_t) i];
        y = pts[3 * (size_t) i + 1];
        z = pts[3 * (size_t) i + 2];
        D[i] = x * x;
        D[n + i] = y * y;
        D[n * 2 + i] = z * z;
        D[n * 3 + i] = 2.0 * y * z;
        D[n * 4 + i] = 2.0 * x * z;
        D[n * 5 + i] = 2.0 * x * y;
        D[n * 6 + i] = 2.0 * x;
        D[n * 7 + i] = 2.0 * y;
        D[n * 8 + i] = 2.0 * z;
        D[n * 9 + i] = 1.0;
    }

    // S = D * D^T (the only use of D)
    Matrix_x_Its_Transpose(S, D, 10, n);
    free(D);
    D = NULL;

    Get_Submatrix(S11, 6, 6, S, 10, 0, 0);
    Get_Submatrix(S12, 6, 4, S, 10, 0, 6);
    Get_Submatrix(S12t, 4, 6, S, 10, 6, 0);
    Get_Submatrix(S22, 4, 4, S, 10, 6, 6);

    for (i = 0; i < 16; i++)
        S22_1[i] = S22[i];
    if (Choleski_LU_Decomposition(S22_1, 4) != 0 || Choleski_LU_Inverse(S22_1, 4) != 0)
        return fit_fail(quiet, "Ellipsoid fit failed (step 1: inverting the 4x4 data sub-matrix): the data is degenerate. "
                               "Rotate the sensor through all directions so the points do not lie on a plane, line or point.\n");

    // Calculate S22a = S22_1 * S12t   4*6 = 4x4 * 4x6   C = AB
    Multiply_Matrices(S22a, S22_1, 4, 4, S12t, 6);

    // Then calculate S22b = S12 * S22a      ( 6x6 = 6x4 * 4x6)
    Multiply_Matrices(S22b, S12, 6, 4, S22a, 6);

    // Calculate SS = S11 - S22b
    for (i = 0; i < 36; i++)
        SS[i] = S11[i] - S22b[i];
    Multiply_Matrices(E, C, 6, 6, SS, 6);

    if (Hessenberg_Form_Elementary(E, SSS, 6) != 0)
        return fit_fail(quiet, "Ellipsoid fit failed (step 2: Hessenberg reduction of the 6x6 matrix): out of memory.\n");

    if (QR_Hessenberg_Matrix(E, SSS, eigen_real, eigen_imag, 6, 100) != 0)
        return fit_fail(quiet, "Ellipsoid fit failed (step 3: eigenvalues of the 6x6 matrix did not converge): "
                               "the data is not usable. Rotate the sensor through all directions and check for bad readings.\n");

    index = 0;
    maxval = eigen_real[0];
    for (i = 1; i < 6; i++) {
        if (eigen_real[i] > maxval) {
            maxval = eigen_real[i];
            index = i;
        }
    }

    v1[0] = SSS[index];
    v1[1] = SSS[index + 6];
    v1[2] = SSS[index + 12];
    v1[3] = SSS[index + 18];
    v1[4] = SSS[index + 24];
    v1[5] = SSS[index + 30];

    // normalize v1
    norm = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] + v1[3] * v1[3] + v1[4] * v1[4] + v1[5] * v1[5]);
    if (!isfinite(norm) || !(norm > 0.0))
        return fit_fail(quiet, "Ellipsoid fit failed (step 4: eigenvector of the 6x6 matrix is zero or not finite): "
                               "the data is degenerate. Rotate the sensor through all directions.\n");
    v1[0] /= norm;
    v1[1] /= norm;
    v1[2] /= norm;
    v1[3] /= norm;
    v1[4] /= norm;
    v1[5] /= norm;

    if (v1[0] < 0.0) {
        v1[0] = -v1[0];
        v1[1] = -v1[1];
        v1[2] = -v1[2];
        v1[3] = -v1[3];
        v1[4] = -v1[4];
        v1[5] = -v1[5];
    }

    // Calculate v2 = S22a * v1      ( 4x1 = 4x6 * 6x1)
    Multiply_Matrices(v2, S22a, 4, 6, v1, 1);

    v[0] = v1[0];
    v[1] = v1[1];
    v[2] = v1[2];
    v[3] = v1[3];
    v[4] = v1[4];
    v[5] = v1[5];
    v[6] = -v2[0];
    v[7] = -v2[1];
    v[8] = -v2[2];
    v[9] = -v2[3];

    Q[0] = v[0];
    Q[1] = v[5];
    Q[2] = v[4];
    Q[3] = v[5];
    Q[4] = v[1];
    Q[5] = v[3];
    Q[6] = v[4];
    Q[7] = v[3];
    Q[8] = v[2];

    U[0] = v[6];
    U[1] = v[7];
    U[2] = v[8];
    for (i = 0; i < 9; i++)
        Q_1[i] = Q[i];
    if (Choleski_LU_Decomposition(Q_1, 3) != 0 || Choleski_LU_Inverse(Q_1, 3) != 0)
        return fit_fail(quiet, "Ellipsoid fit failed (step 5: the fitted quadric is not an ellipsoid, matrix not positive definite): "
                               "the data does not span enough orientations. Rotate the sensor through all directions.\n");

    // Calculate B = Q-1 * U   ( 3x1 = 3x3 * 3x1)
    Multiply_Matrices(B, Q_1, 3, 3, U, 1);
    B[0] = -B[0];     // x-axis combined bias
    B[1] = -B[1];     // y-axis combined bias
    B[2] = -B[2];     // z-axis combined bias

    // First calculate QB = Q * B   ( 3x1 = 3x3 * 3x1)
    Multiply_Matrices(QB, Q, 3, 3, B, 1);

    // Then calculate btqb = BT * QB    ( 1x1 = 1x3 * 3x1)
    Multiply_Matrices(&btqb, B, 1, 3, QB, 1);

    // Calculate hmb = sqrt(btqb - J).
    J = v[9];
    if (!isfinite(btqb) || !isfinite(J) || !(btqb - J > 0.0))
        return fit_fail(quiet, "Ellipsoid fit failed (step 6: field strength term is not positive): "
                               "the data does not describe an ellipsoid. Rotate the sensor through all directions.\n");
    hmb = sqrt(btqb - J);
    if (!isfinite(hmb) || !(hmb > 0.0))
        return fit_fail(quiet, "Ellipsoid fit failed (step 6: field strength is zero or not finite): "
                               "the data does not describe an ellipsoid. Rotate the sensor through all directions.\n");

    // Calculate SQ, the square root of matrix Q
    if (Hessenberg_Form_Elementary(Q, SSSS, 3) != 0)
        return fit_fail(quiet, "Ellipsoid fit failed (step 7: Hessenberg reduction of the 3x3 matrix): out of memory.\n");

    if (QR_Hessenberg_Matrix(Q, SSSS, eigen_real3, eigen_imag3, 3, 100) != 0)
        return fit_fail(quiet, "Ellipsoid fit failed (step 8: eigenvalues of the 3x3 matrix did not converge): "
                               "the data is not usable. Rotate the sensor through all directions and check for bad readings.\n");

    for (i = 0; i < 3; i++) {
        if (!isfinite(eigen_real3[i]) || !isfinite(eigen_imag3[i]) || !(eigen_real3[i] > 0.0)
            || fabs(eigen_imag3[i]) > 1e-9 * fabs(eigen_real3[i]))
            return fit_fail(quiet, "Ellipsoid fit failed (step 9: eigenvalue %d of the 3x3 matrix is not real and positive): "
                                   "the data does not describe an ellipsoid. Rotate the sensor through all directions.\n", i + 1);
    }

    // normalize eigenvectors
    norm1 = sqrt(SSSS[0] * SSSS[0] + SSSS[3] * SSSS[3] + SSSS[6] * SSSS[6]);
    norm2 = sqrt(SSSS[1] * SSSS[1] + SSSS[4] * SSSS[4] + SSSS[7] * SSSS[7]);
    norm3 = sqrt(SSSS[2] * SSSS[2] + SSSS[5] * SSSS[5] + SSSS[8] * SSSS[8]);
    if (!isfinite(norm1) || !isfinite(norm2) || !isfinite(norm3)
        || !(norm1 > 0.0) || !(norm2 > 0.0) || !(norm3 > 0.0))
        return fit_fail(quiet, "Ellipsoid fit failed (step 10: eigenvector of the 3x3 matrix is zero or not finite): "
                               "the data is degenerate. Rotate the sensor through all directions.\n");
    SSSS[0] /= norm1;
    SSSS[3] /= norm1;
    SSSS[6] /= norm1;
    SSSS[1] /= norm2;
    SSSS[4] /= norm2;
    SSSS[7] /= norm2;
    SSSS[2] /= norm3;
    SSSS[5] /= norm3;
    SSSS[8] /= norm3;

    for (i = 0; i < 9; i++)
        Dz[i] = 0.0;
    Dz[0] = sqrt(eigen_real3[0]);
    Dz[4] = sqrt(eigen_real3[1]);
    Dz[8] = sqrt(eigen_real3[2]);

    Multiply_Matrices(vdz, SSSS, 3, 3, Dz, 3);

    Transpose_Square_Matrix(SSSS, 3);

    Multiply_Matrices(SQ, vdz, 3, 3, SSSS, 3);

    // Soft iron / distortion matrix
    for (i = 0; i < 9; i++)
        A_1[i] = SQ[i] * magnitude / hmb;

    // Report nothing unless every value is finite
    for (i = 0; i < 3; i++) {
        if (!isfinite(B[i]))
            return fit_fail(quiet, "Ellipsoid fit failed (step 11: hard-iron offset is not finite): "
                                   "the data is not usable. Rotate the sensor through all directions.\n");
    }
    for (i = 0; i < 9; i++) {
        if (!isfinite(A_1[i]))
            return fit_fail(quiet, "Ellipsoid fit failed (step 11: soft-iron matrix is not finite): "
                                   "the data is not usable. Rotate the sensor through all directions.\n");
    }

    memcpy(B_out, B, sizeof(B));
    memcpy(A_1_out, A_1, sizeof(A_1));
    return 0;
}

int main(int argc, char *argv[]) {
    int nlines = 0;
    double *pts = NULL;
    size_t ptcap = 0;
    char *lbuf = NULL;
    size_t lcap = 0, llen = 0;
    long lineno = 0;
    int rc;
    double B[3], A_1[9];
    int i;

    // Default values
    double magnitude = 1.0;
    char *filename = "../data.txt";
    int reject_outliers = 0;
    int removed = 0;
    int found = 0;
    int json = 0; // output format: 0 = human readable (default), 1 = machine readable (JSON)

    char *strtodptr;
    const char *prog = (argc > 0 && argv[0] != NULL) ? argv[0] : "amc";
    int argi;

    FILE *fp;

    // Portable command-line parser: -m <magnitude>, -f <filename>, -o, -h
    for (argi = 1; argi < argc; argi++) {
        const char *arg = argv[argi];
        if (arg[0] == '-' && arg[1] == 'm') {
            const char *val;
            if (arg[2] != '\0') {
                val = arg + 2; /* attached: -m1.5 */
            } else {
                if (argi + 1 >= argc) {
                    fprintf(stderr, "Option -m requires an argument. " USAGE "\n", prog);
                    return 1;
                }
                val = argv[++argi];
            }
            magnitude = strtod(val, &strtodptr);
            if (strtodptr == val || *strtodptr != '\0' || !isfinite(magnitude) || magnitude <= 0.0) {
                fprintf(stderr, "Invalid magnitude '%s': must be a finite number greater than 0. " USAGE "\n", val, prog);
                return 1;
            }
        } else if (arg[0] == '-' && arg[1] == 'f') {
            if (arg[2] != '\0') {
                filename = (char *) (arg + 2); /* attached: -fdata.txt */
            } else {
                if (argi + 1 >= argc) {
                    fprintf(stderr, "Option -f requires an argument. " USAGE "\n", prog);
                    return 1;
                }
                filename = argv[++argi];
            }
        } else if (strcmp(arg, "-o") == 0) {
            reject_outliers = 1;
        } else if (strcmp(arg, "-j") == 0) {
            json = 1;
        } else if (strcmp(arg, "-h") == 0) {
            printf("\n" USAGE "\n", prog);
            printf("\nOptions:\n");
            printf("  -m <magnitude>   Specify the magnitude (floating-point number)\n");
            printf("  -f <filename>    Specify the filename for processing\n");
            printf("  -o               Reject outlier points (huge or wrong values) and refit\n");
            printf("  -j               Machine readable output (JSON) instead of human readable text\n");
            printf("  -h               Display this help message\n");
            return 0;
        } else if (arg[0] == '-' && arg[1] != '\0') {
            fprintf(stderr, "Invalid option '%s'. " USAGE "\n", arg, prog);
            return 1;
        } else {
            fprintf(stderr, "Unexpected argument '%s'. " USAGE "\n", arg, prog);
            return 1;
        }
    }


    // Print the parsed values
    if (!json) {
        printf("Magnitude: %lf\n", magnitude);
        printf("Filename: %s\n", filename);
    }

    // Open the user-specified file for reading
    fp = fopen(filename, "r");

    // Check if the file was opened successfully
    if (fp == NULL) {
        fprintf(stderr, "Unable to open the file.\n");
        fprintf(stderr, "File opened: %s\n", filename);
        return 1; // Return an error code
    }
// -------------------------------------------------------

    // Single pass: read each line, append valid points to a growable array
    while ((rc = read_line(fp, &lbuf, &lcap, &llen)) == 1) {
        double pv[3];
        const char *q;
        lineno++;
        if (strlen(lbuf) != llen) {
            fprintf(stderr, "%s: line %ld: invalid character in data\n", filename, lineno);
            goto read_fail;
        }
        for (q = lbuf; *q != '\0' && isspace((unsigned char) *q); q++)
            ;
        if (*q == '\0')
            continue; // blank line
        if (!parse_point(lbuf, pv)) {
            fprintf(stderr, "%s: line %ld: expected three finite numbers (x y z) separated by commas, tabs or spaces\n", filename, lineno);
            goto read_fail;
        }
        if ((size_t) nlines >= ptcap) {
            size_t ncap = ptcap ? ptcap * 2 : 256;
            double *np;
            size_t maxcap = (size_t) (INT_MAX / 10);
            if (maxcap > SIZE_MAX / (10 * sizeof(double)))
                maxcap = SIZE_MAX / (10 * sizeof(double));
            if (ncap > maxcap)
                ncap = maxcap;
            if ((size_t) nlines >= ncap) {
                fprintf(stderr, "%s: too many data points\n", filename);
                goto read_fail;
            }
            np = (double *) realloc(pts, ncap * 3 * sizeof(double));
            if (np == NULL) {
                fprintf(stderr, "Out of memory reading %s\n", filename);
                goto read_fail;
            }
            pts = np;
            ptcap = ncap;
        }
        pts[3 * (size_t) nlines] = pv[0];
        pts[3 * (size_t) nlines + 1] = pv[1];
        pts[3 * (size_t) nlines + 2] = pv[2];
        nlines++;
    }
    if (rc < 0) {
        fprintf(stderr, "Out of memory reading %s\n", filename);
        goto read_fail;
    }
    if (ferror(fp)) {
        fprintf(stderr, "Error reading %s\n", filename);
        goto read_fail;
    }
    fclose(fp);
    free(lbuf);
    lbuf = NULL;

    if (nlines < MIN_POINTS) {
        fprintf(stderr, "%s: found %d valid points, at least %d are needed\n", filename, nlines, MIN_POINTS);
        free(pts);
        return 1;
    }

    found = nlines;
    if (!json) {
        printf("Found %d measurements\n", nlines);
        if (!reject_outliers)
            printf("\nOutlier cleanup (-o): off, all %d measurements used. Run with -o to reject outliers and refit.\n", nlines);
    }

    if (reject_outliers) {
        int total = nlines, pass, coarse_removed, fine_passes = 0;

        // Coarse pass first: gross errors would otherwise ruin the first fit.
        if (prefilter_outliers(pts, &nlines) != 0) {
            fprintf(stderr, "Out of memory\n");
            free(pts);
            return 1;
        }
        coarse_removed = total - nlines;
        if (fit_ellipsoid(pts, nlines, magnitude, B, A_1, 0) != 0) {
            free(pts);
            return 1;
        }
        // Fine pass: drop points that do not sit on the fitted sphere, then refit.
        for (pass = 0; pass < MAX_FILTER_PASSES; pass++) {
            double B2[3], A2[9];
            int dropped = residual_filter(pts, &nlines, B, A_1, magnitude);
            if (dropped < 0) {
                fprintf(stderr, "Out of memory\n");
                free(pts);
                return 1;
            }
            if (dropped == 0)
                break;
            fine_passes++;
            if (fit_ellipsoid(pts, nlines, magnitude, B2, A2, 1) != 0) {
                fprintf(stderr, "Outlier filter stopped: refit after removing points failed, keeping the previous fit.\n");
                break;
            }
            memcpy(B, B2, sizeof(B));
            memcpy(A_1, A2, sizeof(A_1));
        }
        removed = total - nlines;
        if (!json) {
            printf("\nOutlier cleanup (-o): on\n");
            printf("  Coarse pass (raw data, gross errors):  %d removed\n", coarse_removed);
            printf("  Fine pass (fit residuals, %d refit%s):  %d removed\n", fine_passes, fine_passes == 1 ? "" : "s", removed - coarse_removed);
            printf("  Total removed: %d\n", removed);
            printf("  Measurements used: %d of %d\n", nlines, total);
        }
    } else if (fit_ellipsoid(pts, nlines, magnitude, B, A_1, 0) != 0) {
        free(pts);
        return 1;
    }
    free(pts);

    if (json) {
        // Machine readable: a single JSON object, nothing else on stdout
        printf("{\"magnitude\":%.17g,\"filename\":", magnitude);
        print_json_string(filename);
        printf(",\"measurements\":%d,\"outliers_rejected\":%d,\"measurements_used\":%d,", found, removed, nlines);
        printf("\"bias\":[%.17g,%.17g,%.17g],\"matrix\":[", B[0], B[1], B[2]);
        for (i = 0; i < 3; i++)
            printf("%s[%.17g,%.17g,%.17g]", i ? "," : "", A_1[i * 3], A_1[i * 3 + 1], A_1[i * 3 + 2]);
        printf("]}\n");
        return 0;
    }

    // Human readable (default)
    // Print hard iron offset
    printf("\nBias (B), subtract from raw readings:\r\n");
    printf("  x: %lf\r\n", B[0]);
    printf("  y: %lf\r\n", B[1]);
    printf("  z: %lf\r\n", B[2]);

    // Print soft iron / distortion matrix
    printf("\nScale and cross-axis matrix (A_1), multiply after subtracting B:\r\n");
    for (i = 0; i < 3; i++)
        printf("  %lf %lf %lf\r\n", A_1[i * 3], A_1[i * 3 + 1], A_1[i * 3 + 2]);

    printf("\nCalibrated = A_1 * (raw - B)\r\n");

    printf("\nTip: run with -j for machine readable output (JSON).\n");

    return 0;

read_fail:
    fclose(fp);
    free(lbuf);
    free(pts);
    return 1;
}
