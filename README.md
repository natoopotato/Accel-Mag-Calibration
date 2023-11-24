# Accelerometer and Magenetometer Calibration


This is a continuation of work by Merlin from sailboatintruments1 who developed the code for this.
Previously a tool called Magneto was distributed however it is no longer available. The math library used mymathlib.com also no longer exists.

With the help of wayback and some other sources I have been able to recover the code and make it available here.

The test data is from Michael Wronas and is available
![here](https://github.com/michaelwro/accelerometer-calibration/tree/main). Also see for more information on implementation.




## Usage

Build and run. You will be prompted for field magnitude and file name.

You can include relative path names to the data files.
Data files should be three column comma delineated with x y z data. See example data from Michael Wronas.

For accelerometer readings measured in G it would be 1 for field magnitude.



## Math.

This technique was originally developed for calibrating Magnetometers, but has been adapted for use with accelerometers.

Magnetic measurements have two types of error hard iron and soft iron distortions.

Hard iron distortions are caused by objects that generate magnetic fields, like speakers or magnetized iron. If such objects are physically attached to the same reference frame as the sensor, they can permanently bias the sensor's output.

On the other hand, soft iron distortions involve alterations in the existing magnetic field, causing it to stretch or distort based on the sensor's orientation relative to the field.

It's worth noting that, in most cases, hard iron distortions have a more substantial impact on the overall measurement error compared to soft iron distortions.


To translate this theory to accelerometer. Hard Iron becomes bias and soft iron is replaced with axis-nonorthoganality.  
Calibration model is based off of the following equation:

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
s_{31} & s_{33} & s_{33} 
\end{bmatrix}
\begin{bmatrix}
a_{x_\mathrm{raw}} -b_x \\
a_{y_\mathrm{raw}} -b_y \\
a_{z_\mathrm{raw}} -b_z 
\end{bmatrix}
```

This model works on the basis of elipsoid curve fitting.
By taking a number of measurements at various angles and orientations, the elipsoid can be fit to the data.
The calibration matrix is then calculated from the elipsoid parameters.
By taking the inverse of this elipsoid it is possible to calibrate the accelerometer data by correcting a given measurement for the elipsoid.

This equation that is then implemented in software is as follows:

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

Where $` x y z `$ are the coordinates of a point on the surface and $` \textbf{Q} `$ is the quadratic form matrix.
This tranlates the general second polynomial of a quadratic into the following matrix form

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

The centre of the ellipsoid can by calculated as the vector

```math
\textbf{b} = -S^{-1}U
```

This is a 1x3 matrix with ellipsoid offsets. (Hard iron compensation)

The ellipsoids semi-axis' a, b, c are give by:

```math
\textbf{a} = \sqrt{\frac{\textbf{b}^T S \textbf{b} - J}{\lambda_1}}
\quad
\textbf{b} = \sqrt{\frac{\textbf{b}^T S \textbf{b} - J}{\lambda_2}}
\quad
\textbf{c} = \sqrt{\frac{\textbf{b}^T S \textbf{b} - J}{\lambda_3}}
```

Where $\lambda_1, \lambda_2, \lambda_3$ are the eigenvalues of the matrix $S$.
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
3. Scale the pricipal axis of the ellipsoid to the unit sphere.
4. Rotate back the sphere to the original orientation of the ellipsoid.
5. Scale the sphere to the unit sphere. i.e. scale to magnitude of 1. Can also be used to scale to a different radius $R$.

For any point $`\mathrm{x}(x, y, z)`$ on the ellipsoid.

```math
\begin{bmatrix}
x_1\\
x_1\\
x_1
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

Where $` \textbf{B} `$ is the bias matrix.

Step 2. The principal axis of the ellipsoid are alligned with the $x,y,z$ axis of the reference frame. This corresponds to applying an inverse rotation by V, that is $V^T$ or $V^{-1}$[^1].
[^1]: This is because the rotation matrix is a square orthogonal matrix, that is $`V^{-1} = V^T`$c

```math
\begin{bmatrix}
x_2\\
x_2\\
x_2
\end{bmatrix}
=
\mathrm{V^T}
\begin{bmatrix}
x\\
y\\
z
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
0 && \sqrt{\frac{a}{b}} && 0\\
0 && 0 && \sqrt{\frac{a}{c}}
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

Step 5. Scale the sphere to the unit sphere. i.e. scale to magnitude of 1.

```math
\begin{bmatrix}
x_5\\
y_5\\
z_5
\end{bmatrix}
=
\frac{1}{a}
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
\frac{1}{a}
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
\mathrm{R}
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
