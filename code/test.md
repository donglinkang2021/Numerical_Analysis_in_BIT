### LU分解

设有线性方程组

$$
Ax = b
$$

其中，$A$为$n \times n$的矩阵，$b$为$n$维列向量，$x$为$n$维列向量。举个具体例子：

$$
\begin{cases}
a_{11}x_1 + a_{12}x_2 +  a_{13}x_3 = b_1 \\
a_{21}x_1 + a_{22}x_2 +  a_{23}x_3 = b_2 \\
a_{31}x_1 + a_{32}x_2 +  a_{33}x_3 = b_3
\end{cases}
$$

其中，$a_{ij}$为矩阵$A$的元素，$b_i$为列向量$b$的元素，$x_i$为列向量$x$的元素。也即：

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2 \\
b_3
\end{bmatrix}
$$

首先，对于第一个方程，我们选取一个常数$l_{11}$遍除第一行，得到：

$$
\frac{a_{11}}{l_{11}}x_1 + \frac{a_{12}}{l_{11}}x_2 +  \frac{a_{13}}{l_{11}}x_3 = \frac{b_1}{l_{11}}
$$

对于上式，我们令

$$
\begin{bmatrix}
\frac{a_{11}}{l_{11}} & \frac{a_{12}}{l_{11}} & \frac{a_{13}}{l_{11}} \\
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
\frac{b_1}{l_{11}} \\
\end{bmatrix}
\\ \Rightarrow
\begin{bmatrix}
u_{11} & u_{12} & u_{13} \\
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
y_1 \\
\end{bmatrix}
$$

其中，$u_{ij}$为矩阵$U$的元素，$y_i$为列向量$y$的元素(待会也会求解它)。

接下来，我们对第二个方程和第三个方程进行第一步的消元，令$l_{21} = \frac{a_{21}}{u_{11}} \quad l_{31} = \frac{a_{31}}{u_{11}}$

$$
\begin{bmatrix}
a_{21} & a_{22} & a_{23} & b_2 \\
a_{31} & a_{32} & a_{33} & b_3
\end{bmatrix}
-
\begin{bmatrix}
l_{21} \\
l_{31} 
\end{bmatrix}
\begin{bmatrix}
u_{11} & u_{12} & u_{13} & y_1 \\
\end{bmatrix}
=
\begin{bmatrix}
0 & a_{22} - l_{21}u_{12} & a_{23} - l_{21}u_{13} & b_2 - l_{21}y_1 \\
0 & a_{32} - l_{31}u_{12} & a_{33} - l_{31}u_{13} & b_3 - l_{31}y_1
\end{bmatrix}
$$

> 注意：
> - 这里的$y_1$是减去上面求解的$y_1$即可
> - 另外上面可以看到一个很有趣的地方，就是l和u的下标，比如$a_{23} - l_{21}u_{13}$把$l_{21}u_{13}$中的1去掉我们就得到了$a_{23}$的下标，这刚好又是我们当前操作元素的位置；

对于等式右边的矩阵（我们重新写成方程的形式，下面对它继续进行消元工作），

$$
\begin{cases}
(a_{22} - l_{21}u_{12})x_2 + (a_{23} - l_{21}u_{13})x_3 = b_2 - l_{21}y_1 \\
(a_{32} - l_{31}u_{12})x_2 + (a_{33} - l_{31}u_{13})x_3 = b_3 - l_{31}y_1
\end{cases}
$$


我们选取第二个常数$l_{22}$遍除第一行，得到：

$$
\begin{bmatrix}
0 & \frac{a_{22} - l_{21}u_{12}}{l_{22}} & \frac{a_{23} - l_{21}u_{13}}{l_{22}}\\
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
\frac{b_2 - l_{21}y_1}{l_{22}} \\
\end{bmatrix}
\\ \Rightarrow
\begin{bmatrix}
0 & u_{22} & u_{23} \\
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
y_2 \\
\end{bmatrix}
$$

继续对第三个方程进行消元，令$l_{32} = \frac{a_{32} - l_{31}u_{12}}{u_{22}}$

$$
\begin{align}
    &\begin{bmatrix}
    0 & a_{32} - l_{31}u_{12} & a_{33} - l_{31}u_{13} & b_3 - l_{31}y_1
    \end{bmatrix}
    -
    \begin{bmatrix}
    l_{32}
    \end{bmatrix}
    \begin{bmatrix}
    0 & u_{22} & u_{23} & y_2
    \end{bmatrix}\\
    =
    &\begin{bmatrix}
    0 & 0 & a_{33} - l_{31}u_{13} - l_{32}u_{23} & b_3 - l_{31}y_1 - l_{32}y_2
    \end{bmatrix}
\end{align}
$$

对于等式右边的矩阵，我们选取第三个常数$l_{33}$遍除第一行，最终得到：

$$
\begin{bmatrix}
0 & 0 & \frac{a_{33} - l_{31}u_{13} - l_{32}u_{23}}{l_{33}}\\
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
\frac{b_3 - l_{31}y_1 - l_{32}y_2}{l_{33}} \\
\end{bmatrix}
\\ \Rightarrow
\begin{bmatrix}
0 & 0 & u_{33} \\
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
y_3 \\
\end{bmatrix}
$$

至此，我们按照消元的一定顺序步骤，得到了$u$和$l$矩阵，也得到了如下的方程组（矩阵形式）：


$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix} 
b_1 \\
b_2 \\
b_3
\end{bmatrix}
\\ \Rightarrow
\begin{bmatrix}
l_{11} & 0 & 0 \\
l_{21} & l_{22} & 0 \\
l_{31} & l_{32} & l_{33}
\end{bmatrix}
\begin{bmatrix}
u_{11} & u_{12} & u_{13} \\
0 & u_{22} & u_{23} \\
0 & 0 & u_{33}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2 \\
b_3
\end{bmatrix}
\\ \Rightarrow
\begin{bmatrix}
l_{11} & 0 & 0 \\
l_{21} & l_{22} & 0 \\
l_{31} & l_{32} & l_{33}
\end{bmatrix}
\begin{bmatrix}
y_1 \\
y_2 \\
y_3
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2 \\
b_3
\end{bmatrix}
\quad
\begin{bmatrix}
u_{11} & u_{12} & u_{13} \\
0 & u_{22} & u_{23} \\
0 & 0 & u_{33}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
y_1 \\
y_2 \\
y_3
\end{bmatrix}
\\
$$

简洁起见，我们统一将$u$矩阵记为$U$，将$l$矩阵记为$L$，则上述解方程的过程可以写为：

$$
Ax = b
\\ \Rightarrow
LUx = b
\\ \Rightarrow
Ly = b \quad \text{（回代法解出y）}
\\ \Rightarrow
Ux = y \quad \text{（回代法解出x）}
$$






