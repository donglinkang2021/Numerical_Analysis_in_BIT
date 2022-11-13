# 解线性方程组的直接法

### 例1

$$
\begin{cases}
3x_1 + 2x_2 + 5x_3 = 6 \\
-x_1 + 4x_2 + 3x_3 = 5 \\
x_1 - x_2 + 3x_3 = 1
\end{cases}
$$

## 普通解法
代码如下：

```python
import numpy as np
A = np.array([[3, 2, 5], [-1, 4, 3], [1, -1, 3]])
b = np.array([6, 5, 1])
x = np.linalg.solve(A, b)
print(x)
```


```python
import numpy as np
A = np.array([[3, 2, 5], [-1, 4, 3], [1, -1, 3]])
b = np.array([6, 5, 1])
x = np.linalg.solve(A, b)
print(x)
```

    [0.5 1.  0.5]


## 高斯若尔当消元法

> 高斯若尔当消元法是一种解线性方程组的方法，它是一种直接法，即不需要迭代求解，而是直接求解。高斯若尔当消元法的基本思想是将线性方程组转化为一个上三角矩阵，然后逐个求解。
> 
> 另注：高斯若尔当消元法是一种直接法，或者说成是一种朴素想法，但如果初学者在网上搜的时候可能会看到一大堆定义，其中不乏这样的说法
> - 高斯若尔当消元法是要消元直到左边矩阵为单位矩阵的方法，[from csdn](https://blog.csdn.net/daduzimama/article/details/120486666)
> - 高斯若尔当消元法是要消元直到左边矩阵为上三角矩阵的方法, [高等代数第三版上册·丘维声·高等教育出版社](https://item.jd.com/12274342.html)
> - 高斯若尔当消元法是写成[A|I]的形式，然后消元直到左边矩阵为单位矩阵的方法, [from baidu baike](https://baike.baidu.com/item/%E9%AB%98%E6%96%AF-%E8%8B%A5%E5%B0%94%E5%BD%93%E6%B6%88%E5%85%83%E6%B3%95/19969775)
> 
> 这里采用最朴素的解法就好（采取说法二），毕竟我们的最终目的只是解方程而已

<!-- 数学公式 -->

$$
\begin{cases}
3x_1 + 2x_2 + 5x_3 = 6 \\
-x_1 + 4x_2 + 3x_3 = 5 \\
x_1 - x_2 + 3x_3 = 1
\end{cases}
$$


对于上面的例子，我们可以将其转化为增广矩阵：

$$
\begin{bmatrix}
3 & 2 & 5 & 6 \\
-1 & 4 & 3 & 5 \\
1 & -1 & 3 & 1
\end{bmatrix}
$$

然后进行初等行变换，得到：

$$
\begin{bmatrix}
3 & 2 & 5 & 6 \\
0 & 14 & 14 & 21 \\
0 & 5 & -4 & 3
\end{bmatrix}
\dots \text{第一步} \\
$$

$$
\begin{bmatrix}
3 & 2 & 5 & 6 \\
0 & 14 & 14 & 21 \\
0 & 0 & 2 & 1
\end{bmatrix}
\dots \text{第二步}\\
$$

最后利用回代法求解，得到：

$$
\begin{cases}
x_3 = \frac{1}{2} \\
x_2 = \frac{21 - 14 * x_3}{14} \\
x_1 = \frac{6 - 2 * x_2 - 5 * x_3}{3} \\
\end{cases}
\Rightarrow
\begin{cases}
x_3 = 0.5 \\
x_2 = 1     \\
x_1 = 0.5   \\
\end{cases}
$$



```python
import numpy as np

# 高斯消去法
def gaussElimination(A, b):
    n = len(b)
    # 前進消去
    for k in range(0, n-1):  # 举例第一列 k = 0
        factor = A[k, k] / A[k+1:, k] # 求出之后每一行的系数
        A[k+1:, :] = A[k, :] - factor[:,np.newaxis] * A[k+1:, :]  #更新之后每一行 np.newaxis的作用是增加一个维度 (2,1)
        b[k+1:] = b[k, :] - factor[:,np.newaxis] * b[k+1:, :]  #更新之后每一行的b

    # 回代
    x = np.zeros(n)
    x[n-1] = b[n-1] / A[n-1, n-1]
    for k in range(n-2, -1, -1):
        x[k] = (b[k] - np.dot(A[k, k+1:], x[k+1:])) / A[k, k]
    return x

A = np.array([[3, 2, 5], [-1, 4, 3], [1, -1, 3]], dtype=float)
b = np.array([[6, 5, 1]], dtype=float)
b = b.T
gaussElimination(A, b)
```




    array([0.5, 1. , 0.5])



随机测试一个方程组

<!-- 数学公式 -->

$$
\begin{cases}
-23x_1 + 11x_2 + x_3 = 0 \\
11x_1 - 3x_2 - 2x_3 = 3 \\
x_1 - 2x_2 + 2x_3 = -1
\end{cases}
$$



```python
A1 = np.array([[-23, 11 , 1], [11, -3, -2], [1, -2, 2]], dtype=float)
b1 = np.array([[0, 3, -1]], dtype=float)
b1 = b1.T
gaussElimination(A1, b1)
```




    array([1., 2., 1.])



### 逐步求解代码如下


```python
A = np.array([[3, 2, 5], [-1, 4, 3], [1, -1, 3]], dtype=float)
b = np.array([[6, 5, 1]], dtype=float)
# 装置b
b = b.T
tmp_f=A[0,0]/A[1:,0] # [-3.  3.]
A[1:,:] = A[0,:] - tmp_f[:,np.newaxis] * A[1:,:] # A = [[ 3.  2.  5.][ 0.  14.  14.][ 0.  5.  -4.]]
b[1:,:] = b[0,:] - tmp_f[:,np.newaxis] * b[1:,:] # b = [[ 6.][ 21.][ 3.]]
```


```python
tmp_f1=A[1,1]/A[2:,1]
A[2:,:] = A[1,:] - tmp_f1[:,np.newaxis] * A[2:,:] # A = [[ 3.  2.  5.][ 0.  14.  14.][ 0.  0.  25.2]]
b[2:,:] = b[1,:] - tmp_f1[:,np.newaxis] * b[2:,:] # b = [[ 6.][ 21.][ 12.6]]
```


```python
# 回代
x = np.zeros((3,1))
x[2,0] = b[2,0]/A[2,2] # x = [[ 0.][ 0.][ 0.5]]
x[1,0] = (b[1,0] - A[1,2]*x[2,0])/A[1,1] # x = [[ 0.][ 1.][ 0.5]]
x[0,0] = (b[0,0] - A[0,1]*x[1,0] - A[0,2]*x[2,0])/A[0,0] # x = [[ 0.5][ 1.][ 0.5]]
print(x)
```

    [[0.5]
     [1. ]
     [0.5]]



```python
print((-4)*(-2.8) + 14)
tmpa = (-4)*(-2.8) + 14
print(3*(-2.8) + 21)
tmpb = 3*(-2.8) + 21
print(tmpb/tmpa)
```

    25.2
    12.600000000000001
    0.5000000000000001


### 所需前置numpy索引的知识


```python
# demo 重新理解numpy的一些操作
import numpy as np
A = np.array([[3, 2, 5], [-1, 4, 3], [1, -1, 3]])
```


```python
# 学习
print(A[0, 1]) # 第一个参数是行，第二个参数是列
print(A[0, 1].shape) # ()
print(A[0, 1].size) # 1
print(A[0, 1].ndim) # 0
print(A[0, 1].dtype) # int32
print(A[0, 1].itemsize) # 4
print(A[0, 1].nbytes) # 4

print(A[0:2, 0:2]) # 第一个参数是行，第二个参数是列，这里是取前两行，前两列 [[3 2] [-1 4]]
print(A[0:2, 0:2].shape) # shape是一个元组，表示矩阵的行数和列数 (2, 2)
print(A[0:2, 0:2].size) # size是一个整数，表示矩阵的元素个数 4
print(A[0:2, 0:2].ndim) # ndim是一个整数，表示矩阵的维数 2
print(A[0:2, 0:2].dtype) # dtype是一个对象，表示矩阵的数据类型 int32
print(A[0:2, 0:2].itemsize) # itemsize是一个整数，表示矩阵中每个元素的字节数 4
print(A[0:2, 0:2].nbytes) # nbytes是一个整数，表示矩阵中所有元素的字节数 16


```

    2
    ()
    0
    1
    int32
    4
    4
    [[ 3  2]
     [-1  4]]
    (2, 2)
    4
    2
    int32
    4
    16


## 例2

$$
\begin{cases}
3x_1 - x_2 + 4x_3 = 7 \\
-x_1 + 2x_2 - 2x_3 = -1 \\
2x_1 - 3x_2 - 2x_3 = 0
\end{cases}
$$

## 普通解法

代码如下：

```python
import numpy as np
A = np.array([[3, -1, 4], [-1, 2, -2], [2, -3, -2]])
b = np.array([7, -1, 0])
x = np.linalg.solve(A, b)
print(x)
```


```python
# 普通解法先确定结果
import numpy as np
A = np.array([[3, -1, 4], [-1, 2, -2], [2, -3, -2]])
b = np.array([7, -1, 0])
x = np.linalg.solve(A, b)
print(x)
```

    [2.  1.  0.5]



```python
# 高斯若儿当消元法验证
import numpy as np
A = np.array([[3, -1, 4], [-1, 2, -2], [2, -3, -2]], dtype=float)
b = np.array([[7, -1, 0]], dtype=float)
b = b.T
gaussElimination(A, b)
```




    array([2. , 1. , 0.5])



## 高斯消元法与克劳特消元法（Gauss and Crout）

<!-- 数学公式 -->

$$
\begin{cases}
3x_1 - x_2 + 4x_3 = 7 \\
-x_1 + 2x_2 - 2x_3 = -1 \\
2x_1 - 3x_2 - 2x_3 = 0
\end{cases}
$$

相比于高斯消元法，克劳特消元法的基本思想其实也差不多，前者是保证L矩阵的对角线元素为1，而后者是保证U矩阵的对角线元素为1。

对于上面的例子，我们可以将其转化为增广矩阵：

$$
\begin{bmatrix}
3 & -1 & 4 & 7 \\
-1 & 2 & -2 & -1 \\
2 & -3 & -2 & 0
\end{bmatrix}
$$


> 对于为什么要进行LU分解，一开始自己也不是很理解，但其实就是将方程组转化为两个方程组，然后分别求解，最后再合并，这样就可以避免矩阵的逆运算，从而提高计算效率。

这里权当复习整理数值分析里的知识，重新写一下LU分解的过程。

### 高斯消元法

我们实际上是将$AX = B$中的$A$分解为$LU$，其中$U$是一个上三角矩阵，$L$是一个下三角矩阵，$L$的对角线元素都是1，

$$
A = LU
$$

$$
\begin{bmatrix}
3 & -1 & 4 \\
-1 & 2 & -2 \\
2 & -3 & -2
\end{bmatrix}
=
\begin{bmatrix}
1 & 0 & 0 \\
-\frac{1}{3} & 1 & 0 \\
\frac{2}{3} & -1.4 & 1
\end{bmatrix}
\begin{bmatrix}
3 & -1 & 4 \\
0 & \frac{5}{3} & -\frac{2}{3} \\
0 & 0 & -5.6
\end{bmatrix}
$$




```python
import numpy as np
A = np.array([[3, -1, 4], [-1, 2, -2], [2, -3, -2]], dtype=float)
b = np.array([[7, -1, 0]], dtype=float)
b = b.T
```


```python
# 高斯消元法
def gauss(A, b,isprint = True):
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    for i in range(n):
        L[i, i] = 1
        for j in range(i, n):
            U[i, j] = A[i, j] - L[i, :i] @ U[:i, j] # @表示矩阵乘法
        for j in range(i, n):
            L[j, i] = (A[j, i] - L[j, :i] @ U[:i, i])/U[i, i]
    y = np.linalg.solve(L, b)
    x = np.linalg.solve(U, y)
    if isprint:
        # 输出解方程的过程
        print("LU分解：")
        print("A = LU")
        print("L = ")
        print(L)
        print("其中我们有：Ly = b，解得")
        print("y = ")
        print(y)
        print("U = ")
        print(U)
        print("其中我们有：Ux = y，解得")
        print("x = ")
        print(x)
    return x

x = gauss(A, b)
```

    LU分解：
    A = LU
    L = 
    [[ 1.          0.          0.        ]
     [-0.33333333  1.          0.        ]
     [ 0.66666667 -1.4         1.        ]]
    其中我们有：Ly = b，解得
    y = 
    [[ 7.        ]
     [ 1.33333333]
     [-2.8       ]]
    U = 
    [[ 3.         -1.          4.        ]
     [ 0.          1.66666667 -0.66666667]
     [ 0.          0.         -5.6       ]]
    其中我们有：Ux = y，解得
    x = 
    [[2. ]
     [1. ]
     [0.5]]


### 克劳特消元法



```python
import numpy as np
A = np.array([[3, -1, 4], [-1, 2, -2], [2, -3, -2]], dtype=float)
b = np.array([[7, -1, 0]], dtype=float)
b = b.T
```


```python
def crout(A, b,isprint = True):
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    for i in range(n):
        L[i, i] = A[i, i] - L[i, :i] @ U[:i, i]
        for j in range(i, n):
            U[i, j] = (A[i, j] - L[i, :i] @ U[:i, j])/L[i, i]
        for j in range(i, n):
            L[j, i] = (A[j, i] - L[j, :i] @ U[:i, i])
    y = np.linalg.solve(L, b)
    x = np.linalg.solve(U, y)
    if isprint:
        # 输出解方程的过程
        print("LU分解：")
        print("A = LU")
        print("L = ")
        print(L)
        print("其中我们有：Ly = b，解得")
        print("y = ")
        print(y)
        print("U = ")
        print(U)
        print("其中我们有：Ux = y，解得")
        print("x = ")
        print(x)
    return x

x = crout(A, b)
```

    LU分解：
    A = LU
    L = 
    [[ 3.          0.          0.        ]
     [-1.          1.66666667  0.        ]
     [ 2.         -2.33333333 -5.6       ]]
    其中我们有：Ly = b，解得
    y = 
    [[2.33333333]
     [0.8       ]
     [0.5       ]]
    U = 
    [[ 1.         -0.33333333  1.33333333]
     [ 0.          1.         -0.4       ]
     [ 0.          0.          1.        ]]
    其中我们有：Ux = y，解得
    x = 
    [[2. ]
     [1. ]
     [0.5]]


### 补充LU分解知识

#### LU分解

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

#### 代码实现


```python
import numpy as np

def lu_decomposition(A):
    """
    LU分解
    :param A: 矩阵
    :return: L, U
    """
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            U[i, j] = A[i, j] - np.dot(L[i, :i], U[:i, j])
        for j in range(i, n):
            if i == j:
                L[i, i] = 1 # 这里的选定将决定是高斯消元法还是克劳特消元法
            else:
                L[j, i] = (A[j, i] - np.dot(L[j, :i], U[:i, i])) / U[i, i]
    return L, U

def forward_substitution(L, b):
    """
    前代法求解Ly = b
    :param L: 下三角矩阵
    :param b: 常数向量
    :return: y
    """
    n = L.shape[0]
    y = np.zeros(n)
    for i in range(n):
        y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]
    return y

def backward_substitution(U, y):
    """
    回代法求解Ux = y
    :param U: 上三角矩阵
    :param y: 常数向量
    :return: x
    """
    n = U.shape[0]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i + 1:], x[i + 1:])) / U[i, i]
    return x

def solve(A, b):
    """
    求解Ax = b
    :param A: 矩阵
    :param b: 常数向量
    :return: x
    """
    L, U = lu_decomposition(A)
    y = forward_substitution(L, b)
    x = backward_substitution(U, y)
    return x

```

#### 测试代码


```python
A = np.array([[3, -1, 4], [-1, 2, -2], [2, -3, -2]], dtype=float)
b = np.array([7, -1, 0], dtype=float) # 注意：这里的b是一维数组，如果要写成矩阵形式即二维向量，需要转置一下
```


```python
lu_decomposition(A)
```




    (array([[ 1.        ,  0.        ,  0.        ],
            [-0.33333333,  1.        ,  0.        ],
            [ 0.66666667, -1.4       ,  1.        ]]),
     array([[ 3.        , -1.        ,  4.        ],
            [ 0.        ,  1.66666667, -0.66666667],
            [ 0.        ,  0.        , -5.6       ]]))




```python
# b = b.T
x = solve(A, b)
print(x)
```

    [2.  1.  0.5]



```python
A = np.array([[3, 2, 5], [-1, 4, 3], [1, -1, 3]])
b = np.array([6, 5, 1])
L,U = lu_decomposition(A)
print(L)
print(U)
y = forward_substitution(L, b)
print(y)
x = backward_substitution(U, y)
print(x)
```

    [[ 1.          0.          0.        ]
     [-0.33333333  1.          0.        ]
     [ 0.33333333 -0.35714286  1.        ]]
    [[3.         2.         5.        ]
     [0.         4.66666667 4.66666667]
     [0.         0.         3.        ]]
    [6.  7.  1.5]
    [0.5 1.  0.5]



```python
-5/14
```




    -0.35714285714285715



#### 小结

##### LU分解

总而言之，LU分解的过程就是将方程组$Ax = b$转化为$LUx = b$的过程，其中$U$是上三角矩阵，$L$是下三角矩阵。

我们逐步完成这个过程的同时也逐渐实现了消元的目的，整个过程简洁如下：

- 自己选定$l_{ii}$这个常数，若 $l_{ii} = 1$那么我们就是用高斯消元法，若$l_{ii} = a_{ii} - \sum_{k=1}^{i-1}l_{ik}u_{ki}$则我们就是用克劳特消元法，

- 其他的常数都是通过消元得到的，如：

$$
\begin{align}
l_{ij} &= \frac{a_{ij} - \sum_{k=1}^{i-1}l_{ik}u_{kj}}{u_{ii}}\\
u_{ij} &= \frac{a_{ij} - \sum_{k=1}^{i-1}l_{ik}u_{kj}}{l_{ii}}\\
\end{align}
$$

LU分解的主要部分写成代码就是这七行：

```python
for i in range(n):
    # L[i, i] = A[i, i] - L[i, :i] @ U[:i, i] 克劳特消元法
    L[i, i] = 1 # 高斯消元法
    for j in range(i, n):
        U[i, j] = (A[i, j] - L[i, :i] @ U[:i, j])/L[i, i]
    for j in range(i, n):
        L[j, i] = (A[j, i] - L[j, :i] @ U[:i, i])/U[i, i]
```

##### 回代法

回代其实就非常简单了

$$
\begin{align}
y_i &= \frac{b_i - \sum_{k=1}^{i-1}l_{ik}y_k}{l_{ii}}\\
x_i &= \frac{y_i - \sum_{k=i+1}^{n}u_{ik}x_k}{u_{ii}}\\
\end{align}
$$

回代法的主要部分写成代码就是这四行：

```python
for i in range(n):
    y[i] = (b[i] - L[i, :i] @ y[:i])/L[i, i]
for i in range(n-1, -1, -1):
    x[i] = (y[i] - U[i, i+1:] @ x[i+1:])/U[i, i]
```

> 注：上面出现过的@符号是矩阵乘法的符号，是python3.5之后才支持的，如果你的python版本低于3.5，可以使用np.dot()代替。

## 例3

$$
\begin{cases}
3x_1 - x_2 + 2x_3 = 7 \\
-x_1 + 2x_2 - 2x_3 = -1 \\
2x_1 - 2x_2 + 4x_3 = 0
\end{cases}
$$

## 普通解法



```python
import numpy as np
A = np.array([
    [3, -1, 2], 
    [-1, 2, -2], 
    [2, -2, 4]
    ])
b = np.array([
    7, -1, 0
    ])
x = np.linalg.solve(A, b)
print(x)
```

    [ 3.5  -1.   -2.25]


## 平方根法(Cholesky decomposition)

当A是一个对称正定矩阵，我们在选取$l_{ii}$的值时就可以做一些对称的操作，比如直接令$l_{ii} = u_{ii}$，这样就可以减少一些计算量。

在平方根法中，我们选取$l_{ii} = \sqrt{a_{ii} - \sum_{k=1}^{i-1}l_{ik}^2}$作为$l_{ii}$，这样我们就可以将$A$分解为$LL^T$，其中$U = L^T$。

$$
\begin{align}
A &= LL^T\\
A &= \begin{bmatrix}
l_{11} & 0 & 0 \\
l_{21} & l_{22} & 0 \\
l_{31} & l_{32} & l_{33}
\end{bmatrix}
\begin{bmatrix}
l_{11} & l_{21} & l_{31} \\
0 & l_{22} & l_{32} \\
0 & 0 & l_{33}
\end{bmatrix}\\
&= \begin{bmatrix}
l_{11}^2 & l_{11}l_{21} & l_{11}l_{31} \\
l_{21}l_{11} & l_{21}^2 + l_{22}^2 & l_{21}l_{31} + l_{22}l_{32} \\
l_{31}l_{11} & l_{31}l_{21} + l_{32}l_{22} & l_{31}^2 + l_{32}^2 + l_{33}^2
\end{bmatrix}\\
\end{align}
$$

### 代码实现



```python
def cholesky(A):
    n = A.shape[0]
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] = np.sqrt(A[i, i] - L[i, :i] @ L[i, :i])
        for j in range(i+1, n):
            L[j, i] = (A[j, i] - L[j, :i] @ L[i, :i])/L[i, i]
    return L
```

### 测试代码


```python
A = np.array([
    [3, -1, 2], 
    [-1, 2, -2], 
    [2, -2, 4]
    ])
cholesky(A)
```




    array([[ 1.73205081,  0.        ,  0.        ],
           [-0.57735027,  1.29099445,  0.        ],
           [ 1.15470054, -1.03279556,  1.26491106]])




```python
b = np.array([
    7, -1, 0
    ])
L = cholesky(A)
y = np.linalg.solve(L, b)
print(y)
x = np.linalg.solve(L.T, y)
print(x)
```

    [ 4.04145188  1.03279556 -2.84604989]
    [ 3.5  -1.   -2.25]


### 补充numpy中矩阵的一系列自带函数


```python
import numpy as np
A = np.array([
    [3, -1, 2], 
    [-1, 2, -2], 
    [2, -2, 4]
    ])
# 求A的转置
print(A.T)
'''
[[ 3 -1  2]
 [-1  2 -2]
 [ 2 -2  4]]
 '''
# 求A的逆矩阵 
A_inv = np.linalg.inv(A)
print(A_inv)
'''
[[ 0.5    0.    -0.25 ]
 [ 0.     1.     0.5  ]
 [-0.25   0.5    0.625]]
'''
# 求A的行列式
A_det = np.linalg.det(A)
print(A_det)
'''
8.000000000000002
'''
# 求A的特征值和特征向量
A_eig = np.linalg.eig(A)
print(A_eig)
'''
(array([6.61185871, 1.65867531, 0.72946598]), 
array([[ 0.52548211, -0.83393758, -0.16857242],
       [-0.43184431, -0.43214537,  0.7916823 ],
       [ 0.73306142,  0.34321785,  0.58721586]]))
'''
# 求A的秩
A_rank = np.linalg.matrix_rank(A)
print(A_rank) # 3
# 求A的范数
A_norm = np.linalg.norm(A)
print(A_norm) # 6.855654600401044
# 求A的条件数
A_cond = np.linalg.cond(A)
print(A_cond) # 9.063971310181588
# 求A的奇异值分解
A_svd = np.linalg.svd(A)
print(A_svd)
'''
(array([
    [-0.52548211, -0.83393758, -0.16857242],
    [ 0.43184431, -0.43214537,  0.7916823 ],
    [-0.73306142,  0.34321785,  0.58721586]
    ]), 
 array([6.61185871, 1.65867531, 0.72946598]), 
 array([
    [-0.52548211,  0.43184431, -0.73306142],
    [-0.83393758, -0.43214537,  0.34321785],
    [-0.16857242,  0.7916823 ,  0.58721586]
    ]))
'''
# 求A的QR分解
A_qr = np.linalg.qr(A)
print(A_qr)
'''
(
array([[-0.80178373, -0.5179324 , -0.2981424 ],
       [ 0.26726124, -0.75697812,  0.59628479],
       [-0.53452248,  0.39840954,  0.74535599]]), 
array([[-3.74165739,  2.40535118, -4.27617987],
       [ 0.        , -1.79284291,  2.07172959],
       [ 0.        ,  0.        ,  1.19256959]])
)
'''
```

## 总结

在数值分析这门课中的解线性方程组的直接法，个人认为掌握LU分解一个点就够了；平方根法的话，如果你的矩阵是对称正定的，那么可以考虑使用平方根法，因为它的计算量比LU分解小很多。

此外，我们还学习了其它的解方程的直接法，比如列主元素法，全主元素法，个人觉得这些方法都是在高斯消元法的基础上做了一些优化，所以我就不再赘述了。（因为优化程度不太大，只是选取一个最优的参数来给$l_{ii}$操作而已）

还有针对特殊的矩阵，比如对称矩阵，我们还可以使用平方根法，这个方法在上面已经有所介绍了。针对三对角矩阵，我们可以使用Thomas算法（课本写作追赶法），但还是那句话，都是在高斯消元法的基础上做了一些优化而已；

下一章我们将学习迭代法来求解线性方程组。

### 一些测试


```python
A = np.array([
    [2, -1, 0, 0],
    [-1, 2, -1, 0],
    [0, -1, 2, -1],
    [0, 0, -1, -2]
    ])
b = np.array([
    0, 1, 0, 2.5
    ])
# 求解线性方程组
x = np.linalg.solve(A, b)
print(x) 
```

    [ 0.22727273  0.45454545 -0.31818182 -1.09090909]



```python
lu_decomposition(A)
```




    (array([[ 1.        ,  0.        ,  0.        ,  0.        ],
            [-0.5       ,  1.        ,  0.        ,  0.        ],
            [ 0.        , -0.66666667,  1.        ,  0.        ],
            [ 0.        ,  0.        , -0.75      ,  1.        ]]),
     array([[ 2.        , -1.        ,  0.        ,  0.        ],
            [ 0.        ,  1.5       , -1.        ,  0.        ],
            [ 0.        ,  0.        ,  1.33333333, -1.        ],
            [ 0.        ,  0.        ,  0.        , -2.75      ]]))




```python
tmp = np.array([
    [2, -1, 0, 0],
    [0,1.5, -1, 0],
    [0, 0, 4/3, -1],
    [0, 0, 0, -2.75]
    ])
tmpy = np.array([
    0, 1, 2/3, 3
    ])
np.linalg.solve(tmp, tmpy)
```




    array([ 0.22727273,  0.45454545, -0.31818182, -1.09090909])

