# Accelerometer and Magnetometer Calibration


This is a continuation of work by Merlin from sailboatinstruments1, who wrote the original code.
The code came from a tool called Magneto. It used to be distributed but isn't available anymore. The math library it used, mymathlib.com, is gone too.

I recovered the code with the help of Wayback and some other sources and put it here.

The test data is from Michael Wronas and is available
[here](https://github.com/michaelwro/accelerometer-calibration/tree/main). His repo also has more on the implementation.

He also has a good [video](https://youtu.be/-1tmYPE7MAQ) about this.


## Usage

Build with CMake, then run the `accel_calc` binary:

```
cmake -S . -B build
cmake --build build
./build/accel_calc -m 1 -f data.txt
./build/accel_calc -h
```

With a multi-config generator such as MSVC the binary lands in `build/Debug` or `build/Release` instead.

- `-m <magnitude>` is the target field magnitude, in the units you want the output in. It must be a finite number greater than 0, and defaults to 1.0.
  For an accelerometer use 1 (readings in g) or 9.80665 (m/s^2).
  For a magnetometer use the local Earth-field magnitude (for example in µT), or 1 for a unit sphere.
  Both `-m 9.81` and `-m9.81` work.
- `-f <filename>` is the data file to read. Both `-f data.txt` and `-fdata.txt` work, and a relative path is resolved against the current working directory.
  The default is `../data.txt`, which works when you run the program from `build/`.
- `-o` rejects outlier points and refits (see Outlier rejection below). Off by default.
- `-j` switches to machine readable output (a single JSON object). The default is human readable text.
- `-h` prints the help text and exits.

### Input format

One point per line: exactly three finite numbers x y z, separated by spaces or tabs, with at most one comma between numbers.
Blank lines are ignored and CRLF line endings are fine.
Header rows, comments, extra columns and `nan`/`inf` values stop the run with a `line N: ...` error.
The program needs at least 10 points, but a good fit needs far more (see Collecting data below). See the example data from Michael Wronas.

### Output

Everything goes to stdout, in this order:

1. `Magnitude:`, `Filename:` and `Found N measurements` lines, then an "Outlier cleanup" block. With `-o` it reports how many points the coarse pass and the fine pass removed, the number of refits, the total removed and how many measurements were used. Without `-o` it says cleanup is off and how to turn it on.
2. A labelled "Bias (B)" block with three lines, `x:`, `y:` and `z:`.
3. A labelled "Scale and cross-axis matrix (A_1)" block with the 3x3 matrix, three values per row, then a reminder of the formula.

B and $A_1$ are printed with 6 decimals and `\r\n` line endings, followed by a tip line about `-j`. Errors go to stderr and the exit status is 1.

With `-j` stdout holds only one JSON object, with no header lines and no tip:

```
{"magnitude":1,"filename":"data.txt","measurements":500,"outliers_rejected":0,"measurements_used":500,"bias":[bx,by,bz],"matrix":[[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]]}
```

Numbers are printed with 17 significant digits. Errors still go to stderr with exit status 1.

Apply the result to a raw reading as

```math
\mathrm{calibrated} = A_1 \, (\mathrm{raw} - B)
```

B is in raw units and $A_1$ is in target units per raw unit. With raw counts $A_1$ can be very small and 6 decimals loses precision, so use scaled units (g, m/s^2, µT) in the data file.

### Collecting data

Use hundreds of points that cover all orientations, not just the minimum of 10.
For an accelerometer, hold it still in each pose so that only gravity is measured.
For a magnetometer, keep it away from nearby iron and magnets and sweep it through a full sphere of orientations.

### Outlier rejection

One huge or wrong reading can wreck the fit, because the fit squares every coordinate. With `-o` the program filters in two passes, both using the median and median absolute deviation (MAD) instead of the mean, since a mean is pulled around by the bad points:

1. Coarse pass on the raw data. It finds the centre as the per-axis median and drops points whose distance from that centre is more than 5 robust standard deviations above the median distance. This removes gross errors such as `-3000 4000 1e5` before the first fit.
2. Fine pass on the fit. It calibrates every point, measures how far its length is from the target magnitude `-m`, and drops points more than 3.5 robust standard deviations from the median error (never tighter than 0.5% of the magnitude). It then refits, up to 5 times, until no more points are dropped.

The filter never drops below 10 points. If it would, it stops with a warning on stderr and keeps what it has. It also removes ordinary noisy readings that sit far from the sphere, so on clean data like `data.txt` it can reject a few points.

It cannot rescue data where most points are bad or the orientations are not covered. The first fit is then already wrong and the residuals mislead.

### Errors

Each of these prints a message to stderr and exits with status 1.

| Message | Meaning | Fix |
| --- | --- | --- |
| `Invalid option ...`, `Unexpected argument ...` | Unknown option or stray positional argument. | Use only `-m`, `-f` and `-h`. |
| `Option -m requires an argument`, `Option -f requires an argument` | The option is last on the command line with no value. | Add the value. |
| `Invalid magnitude ...` | `-m` is not a number, or is not finite and greater than 0. | Pass a positive number. |
| `Unable to open the file.` | The file does not exist or cannot be read. | Check the path against the current directory. |
| `<file>: line N: expected three finite numbers ...` | Line N is not exactly three finite numbers (header, comment, extra column, nan/inf). | Fix or remove that line. |
| `<file>: line N: invalid character in data` | Line N contains a NUL byte. | Check the file is plain text. |
| `<file>: found N valid points, at least 10 are needed` | Too few points. | Collect more data. |
| `Error reading ...`, `Out of memory ...`, `too many data points` | I/O failure, or the file is too large. | Retry, or use a smaller file. |
| `Ellipsoid fit failed (step 1 ...)` | The 4x4 data sub-matrix cannot be inverted: the points lie on a plane, line or point. | Rotate the sensor through all directions. |
| `Ellipsoid fit failed (step 2 ...)` or `(step 7 ...)` | Hessenberg reduction ran out of memory. | Free memory and retry. |
| `Ellipsoid fit failed (step 3 ...)` or `(step 8 ...)` | The eigenvalue iteration did not converge. | Check for bad readings and cover all orientations. |
| `Ellipsoid fit failed (step 4 ...)` or `(step 10 ...)` | An eigenvector is zero or not finite: degenerate data. | Cover all orientations. |
| `Ellipsoid fit failed (step 5 ...)` | The fitted quadric is not an ellipsoid (3x3 block not positive definite). | Cover more orientations and remove outliers. |
| `Ellipsoid fit failed (step 6 ...)` | The field-strength term `b^T S b - J` is not positive or finite, so the data does not describe an ellipsoid. | Cover all orientations and remove outliers. |
| `Ellipsoid fit failed (step 9 ...)` | An eigenvalue of the 3x3 block is not real and positive. | Cover all orientations and remove outliers. |
| `Ellipsoid fit failed (step 11 ...)` | The result B or $A_1$ is not finite. | Check the data for bad readings. |

### Limitations

The fitted matrix is symmetric, so it corrects offset, scale and cross-axis skew, but it cannot recover a pure rotation between the sensor frame and the body frame.
The fit minimizes algebraic distance rather than geometric distance, so outliers and poorly distributed data hurt it. Use `-o` to reject outliers.
It uses the ellipsoid-specific constraint from Li & Griffiths, which forces the result to be an ellipsoid and not some other quadric.



## The Math

This technique was originally developed for calibrating Magnetometers, but has been adapted for use with accelerometers.

**Notation.** A few letters are reused below:

- $\textbf{S}$ in the calibration model is the calibration matrix. Further down, $S$ is the 3x3 quadric block $`\begin{bmatrix} A & D & E\\ D & B & F\\ E & F & C \end{bmatrix}`$ of the ellipsoid. The code's `Q` is this block, and the code's $A_1$ is the matrix applied to the data, $\frac{R}{\sqrt{\textbf{b}^T S \textbf{b} - J}}\sqrt{S}$.
- $\textbf{b}$ (bold) is the ellipsoid centre, the hard-iron bias vector. It is the code's `B`. The plain $a, b, c$ are the semi-axis lengths, and $\textbf{B}$ in step 1 is the same vector as $\textbf{b}$.
- $\sqrt{\textbf{b}^T S \textbf{b} - J}$ is the code's `hmb`.
- $R$ is the target magnitude, given by `-m`.

Magnetic measurements have two types of error: hard iron and soft iron distortion.

Hard iron distortion comes from objects that generate magnetic fields, like speakers or magnetized iron. If such an object is attached to the same reference frame as the sensor, it permanently biases the sensor's output.

Soft iron distortion alters the existing magnetic field, stretching or warping it depending on the sensor's orientation relative to the field.

In most cases hard iron distortion adds more error than soft iron distortion.


For an accelerometer, hard iron becomes bias and soft iron becomes axis non-orthogonality.
The calibration model is based on the following equation:

```math
\begin{bmatrix}
a_x \\
a_y\\
a_z
\end{bmatrix}_\mathrm{calibrated}
\quad =
\begin{bmatrix}
s_{11} & s_{12} & s_{13} \\
s_{21} & s_{22} & s_{23} \\
s_{31} & s_{32} & s_{33} 
\end{bmatrix}
\begin{bmatrix}
a_{x_\mathrm{raw}} -b_x \\
a_{y_\mathrm{raw}} -b_y \\
a_{z_\mathrm{raw}} -b_z 
\end{bmatrix}
```

The model fits an ellipsoid to measurements taken at many angles and orientations.
The calibration matrix comes from the ellipsoid parameters, and inverting the ellipsoid corrects a raw measurement.

The software implements this equation:

```math
\begin{bmatrix}
a_x \\
a_y\\
a_z
\end{bmatrix}_{\mathrm{calibrated}}
\quad
=
\textbf{S}^{-1}
\begin{bmatrix}
a_{x_\mathrm{raw}} -b_x \\
a_{y_\mathrm{raw}} -b_y \\
a_{z_\mathrm{raw}} -b_z 
\end{bmatrix}
```

Where $\textbf{S}$ is the calibration matrix and $\textbf{b}$ is the offset vector.


Full general equation of an ellipsoid:

```math
Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
```

Using the matrix form of the quadratic surface.

```math
\textbf{x}
=
\begin{bmatrix}
x\\
y\\
z\\
1
\end{bmatrix}
\quad
\textbf{Q}
=
\begin{bmatrix}
A & D & E & G\\
D & B & F & H\\
E & F & C & I\\
G & H & I & J
\end{bmatrix}
```

Where $x y z$ are the coordinates of a point on the surface and $\textbf{Q}$ is the quadratic form matrix.
This translates the general second-degree polynomial into the following matrix form

```math
\textbf{x}^T \textbf{Q} \textbf{x} = 0
```

By defining the following matrices

```math
S =
\begin{bmatrix}
A & D & E\\
D & B & F\\
E & F & C
\end{bmatrix}
\quad
U =
\begin{bmatrix}
G\\
H\\
I
\end{bmatrix}
```

The centre of the ellipsoid can be calculated as the vector

```math
\textbf{b} = -S^{-1}U
```

This is a 3x1 column vector with the ellipsoid offsets. (Hard iron compensation)

The ellipsoid's semi-axes a, b, c are given by:

```math
\textbf{a} = \sqrt{\frac{\textbf{b}^T S \textbf{b} - J}{\lambda_1}}
\quad
\textbf{b} = \sqrt{\frac{\textbf{b}^T S \textbf{b} - J}{\lambda_2}}
\quad
\textbf{c} = \sqrt{\frac{\textbf{b}^T S \textbf{b} - J}{\lambda_3}}
```

Where $`\lambda_1, \lambda_2, \lambda_3`$ are the eigenvalues of the matrix $S$.
The corresponding normalised eigenvectors are the semi-axis directions.
The semi-axis directions are the columns of the matrix $V$.

```math
V =
\begin{bmatrix}
\mathrm{v}_{1_{\mathrm{3x1}}} & \mathrm{v}_{2_{\mathrm{3x1}}} & \mathrm{v}_{3_{\mathrm{3x1}}}
\end{bmatrix}
```

#### Procedure to map a general ellipsoid to a unit sphere centered at the origin

1. Translate the ellipsoid to the origin.
2. Align the principal axis of the ellipsoid with the coordinate axis.
3. Scale the principal axis of the ellipsoid to the unit sphere.
4. Rotate back the sphere to the original orientation of the ellipsoid.
5. Scale the sphere to the target radius $R$ (the magnitude given with `-m`). $R = 1$ gives the unit sphere.

For any point $\mathrm{x}(x, y, z)$ on the ellipsoid.

```math
\begin{bmatrix}
x_1\\
y_1\\
z_1
\end{bmatrix}
=
\begin{bmatrix}
x\\
y\\
z
\end{bmatrix}
-
\textbf{B}
=
\begin{bmatrix}
x - b_x\\
y - b_y\\
z - b_z
\end{bmatrix}
```

Where $\textbf{B}$ is the bias vector.

Step 2. The principal axes of the ellipsoid are aligned with the $x,y,z$ axes of the reference frame. This corresponds to applying an inverse rotation by V, that is $V^T$ or $V^{-1}$[^1].
[^1]: $V^T = V^{-1}$ holds because $S$ is symmetric, so its eigenvectors are orthonormal and $V$ is an orthogonal matrix (a rotation, up to the sign of a column).

```math
\begin{bmatrix}
x_2\\
y_2\\
z_2
\end{bmatrix}
=
\mathrm{V^T}
\begin{bmatrix}
x_1\\
y_1\\
z_1
\end{bmatrix}
```

Step 3. The semi-axis a, b, c are scaled to the length of the semi-axis a, transforming the ellipsoid into a sphere of radius a.

```math
\begin{bmatrix}
x_3\\
y_3\\
z_3
\end{bmatrix}
=
\begin{bmatrix}
1 && 0 && 0\\
0 && \frac{a}{b} && 0\\
0 && 0 && \frac{a}{c}
\end{bmatrix}
\begin{bmatrix}
x_2\\
y_2\\
z_2
\end{bmatrix}
=
\begin{bmatrix}
1 && 0 && 0\\
0 && \sqrt{\frac{\lambda_2}{\lambda_1}} && 0\\
0 && 0 && \sqrt{\frac{\lambda_3}{\lambda_1}}
\end{bmatrix}
\begin{bmatrix}
x_2\\
y_2\\
z_2
\end{bmatrix}
=
\frac{1}{\sqrt{\lambda_1}}
\begin{bmatrix}
\sqrt{\lambda_1} && 0 && 0\\
0 && \sqrt{\lambda_2} && 0\\
0 && 0 && \sqrt{\lambda_3}
\end{bmatrix}
\begin{bmatrix}
x_2\\
y_2\\
z_2
\end{bmatrix}
```

Step 4. Rotate the sphere back to the original orientation of the ellipsoid. This corresponds to applying the rotation matrix V.

```math
\begin{bmatrix}
x_4\\
y_4\\
z_4
\end{bmatrix}
=
\mathrm{V}
\begin{bmatrix}
x_3\\
y_3\\
z_3
\end{bmatrix}
```

Step 5. Scale the sphere of radius $a$ to the target radius $R$ ($R = 1$ gives the unit sphere).

```math
\begin{bmatrix}
x_5\\
y_5\\
z_5
\end{bmatrix}
=
\frac{R}{a}
\begin{bmatrix}
x_4\\
y_4\\
z_4
\end{bmatrix}
```

The final equation is then in various forms:

```math
\begin{bmatrix}
x_5\\
y_5\\
z_5
\end{bmatrix}
=
\frac{R}{a\sqrt{\lambda_1}}
\mathrm{V}
\begin{bmatrix}
\sqrt{\lambda_1} && 0 && 0\\
0 && \sqrt{\lambda_2} && 0\\
0 && 0 && \sqrt{\lambda_3}\\
\end{bmatrix}
\mathrm{V^T}
\begin{bmatrix}
x - b_x\\
y - b_y\\
z - b_z
\end{bmatrix}
```

```math
\begin{bmatrix}
x_5\\
y_5\\
z_5
\end{bmatrix}
=
R
\sqrt{\frac{\lambda_1}{\textbf{b}^T S \textbf{b} - J}}
\sqrt{\frac{1}{\lambda_1}}\mathrm{V}
\begin{bmatrix}
\sqrt{\lambda_1} && 0 && 0\\
0 && \sqrt{\lambda_2} && 0\\
0 && 0 && \sqrt{\lambda_3}
\end{bmatrix}
\mathrm{V^T}
\begin{bmatrix}
x - b_x\\
y - b_y\\
z - b_z
\end{bmatrix}
```

```math
\begin{bmatrix}
x_5\\
y_5\\
z_5
\end{bmatrix}
=
\frac{R}{\sqrt{\textbf{b}^T S \textbf{b} - J}}
\mathrm{V}
\begin{bmatrix}
\sqrt{\lambda_1} && 0 && 0\\
0 && \sqrt{\lambda_2} && 0\\
0 && 0 && \sqrt{\lambda_3}
\end{bmatrix}
\mathrm{V^T}
\begin{bmatrix}
x - b_x\\
y - b_y\\
z - b_z
\end{bmatrix}
```

```math
\begin{bmatrix}
x_5\\
y_5\\
z_5
\end{bmatrix}
=
\frac{R}{\sqrt{\textbf{b}^T S \textbf{b} - J}}
\sqrt{\textbf{S}}
\begin{bmatrix}
x - b_x\\
y - b_y\\
z - b_z
\end{bmatrix}
```

## References

- Q. Li and J. G. Griffiths, "Least Squares Ellipsoid Specific Fitting", Proc. Geometric Modeling and Processing (GMP), 2004.
- Magneto v1.2 by Merlin (sailboatinstruments1), whose code this continues.
- The matrix routines in `mathlib/`, from mymathlib.com (the source files do not name their origin in their headers).
