# 解线性方程组的迭代法

## 前言

> - 本文主要用两个简单的例子来介绍了解线性方程组的三种迭代法的原理和实现方法：第一个例子供我们去学习，而第二个例子供我们去验证。
> - 还另外介绍了一些范数的知识，以及如何用Python来计算范数。
> - 另：本文的代码实现全部基于Python。

### 例1

求解方程组

$$
\begin{cases}
20x_1 + 2x_2 + 3x_3 = 24 \\
x_1 + 8x_2 + x_3 = 12 \\
2x_1 - 3x_2 + 15x_3 = 30
\end{cases}
$$

写成矩阵形式

$$
\begin{bmatrix}
20 & 2 & 3 \\
1 & 8 & 1 \\
2 & -3 & 15
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
24 \\
12 \\
30
\end{bmatrix}
\\
令A =
\begin{bmatrix}
20 & 2 & 3 \\
1 & 8 & 1 \\
2 & -3 & 15
\end{bmatrix}
,
X =
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
,
B =
\begin{bmatrix}
24 \\
12 \\
30
\end{bmatrix}
\\ 则有AX = B
$$



```python
import numpy as np
import pandas as pd

A = np.array([
    [20, 2, 3], 
    [1, 8, 1], 
    [2, -3, 15]
    ])
B = np.array([
    [24], 
    [12], 
    [30]
    ])
X = np.array([
    [0], 
    [0], 
    [0]
    ])
```


```python
np.linalg.solve(A, B)
```




    array([[0.76735381],
           [1.13840976],
           [2.12536811]])



### 朴素想法

建立迭代公式

$$
AX = B
\\ \Rightarrow
0 = B - AX
\\ \Rightarrow
X = (I-A)X + B
$$



```python
def simple_iteration(A, B ,n = 10):
    X = np.array([[0], [0], [0]])
    # 单位矩阵
    I = np.eye(A.shape[0])
    for i in range(n):
        X = np.dot(I-A, X) + B
    return X

simple_iteration(A, B) # 10次迭代发现下面结果不收敛,果断舍弃,下面我们将介绍三种收敛的迭代方法
```




    array([[-1.70908945e+13],
           [-1.69858915e+12],
           [-4.97527196e+12]])



## 简单迭代法

### 雅可比迭代法

#### 数学推理

对于一个如下的n阶线性方程组

$$
\begin{cases}
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1 \\
a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2 \\
 \vdots \\
a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nn}x_n = b_n
\end{cases}
$$

写成如下求和形式

$$
\sum_{j=1}^n a_{ij}x_j = b_i, i=1,2,\cdots,n \\
$$

则对于每一列来说，都有

$$
x_i = (b_i - \sum_{k=1,k\neq i}^n a_{ij}x_k)/a_{ii}, i=1,2,\cdots,n
$$

写成矩阵形式

$$
X = (B + (D-A)X)D^{-1} \\
其中D为对角矩阵，D_{ii} = a_{ii}, i=1,2,\cdots,n \\
$$

令$G = D^{-1}(D-A) = I - D^{-1}A, F=D^{-1}B$，则有

$$
X = GX + F
$$

这就是我们的雅可比迭代公式。（又叫简单迭代法）



#### 代码实现


```python
import numpy as np

def jacobi(A, B, X0, eps=1e-6, maxIter=10):
    n = len(A)
    X = X0
    D = np.diag(A) # 对角线元素, 返回一维数组
    R = np.diagflat(D) - A # 对角线元素为0的矩阵, diagflat()函数用于将一维数组转换为对角矩阵
    for i in range(maxIter):
        print('第{}次迭代'.format(i))
        print('X = {}'.format(X))
        X = (B + np.dot(R, X0)) / D
        if np.linalg.norm(X - X0) < eps:
            print('迭代次数：', i)
            break
        X0 = X
    return X


A = np.array([[20, 2, 3], [1, 8, 1], [2, -3, 15]])
B = np.array([24, 12, 30])
X0 = np.array([1, 1, 2])
X = jacobi(A, B, X0)
print(X)
```

    第0次迭代
    X = [1 1 2]
    第1次迭代
    X = [0.8        1.125      2.06666667]
    第2次迭代
    X = [0.7775     1.14166667 2.11833333]
    第3次迭代
    X = [0.76808333 1.13802083 2.12466667]
    第4次迭代
    X = [0.76749792 1.13840625 2.12519306]
    第5次迭代
    X = [0.76738042 1.13841363 2.12534819]
    第6次迭代
    X = [0.76735641 1.13840892 2.12536534]
    第7次迭代
    X = [0.76735431 1.13840978 2.1253676 ]
    迭代次数： 7
    [0.76735388 1.13840976 2.12536805]
    


```python
A @ X # 检验结果
```




    array([24.00000054, 12.00000009, 29.9999996 ])



#### 补充知识——范数

##### 向量范数

范数是向量空间中的一个函数，它将向量映射到实数域上。范数的定义是：对于向量空间中的任意两个向量，都有一个非负的实数与之对应，这个实数称为这两个向量的范数。范数的性质：

1. 非负性：$||x|| \geq 0$
2. 对任意实数$\alpha$，有$||\alpha x|| = |\alpha||x||$
3. 对任意向量$X,Y \in \mathbb{R}^n$，有$||X+Y|| \leq ||X|| + ||Y||$

一些常用范数：

1. $L_1$范数：$||x||_1 = \sum_{i=1}^n |x_i|$
2. $L_2$范数：$||x||_2 = \sqrt{\sum_{i=1}^n x_i^2}$
3. $L_\infty$范数：$||x||_\infty = \max_{i=1,2,\cdots,n}|x_i|$

其中$L_p$范数定义为

$$
||x||_p = \left(\sum_{i=1}^n |x_i|^p\right)^{1/p}
$$

当不需要区分范数时，我们一般用$||x||$表示$L_2$范数。（书上说的是泛指任何一种向量范数，但是我觉得$L_2$更常用）

##### 矩阵范数

常用的矩阵范数有：

1. $L_1$范数：$||A||_1 = \max_{1\leq j\leq n}\sum_{i=1}^n |a_{ij}|$ （列和范数）
2. $L_\infty$范数：$||A||_\infty = \max_{1\leq i\leq n}\sum_{j=1}^n |a_{ij}|$ （行和范数）
3. $Frobenius$范数：$||A||_F = \sqrt{\sum_{i=1}^n\sum_{j=1}^n a_{ij}^2}$ （矩阵元素平方和的平方根）
4. $L_2$范数：$||A||_2 = \max_{\|x\|=1}||Ax||_2$

其中$L_p$范数定义为

$$
||A||_p = \max_{\|x\|=1}||Ax||_p
$$

其中$\|x\|$表示$x$的$L_2$范数，x是一个向量。

书本上的二范数定义为

$$
||A||_2 = \sqrt{\lambda_{max}(A^*A)} \\
其中\lambda_{max}是矩阵A的最大特征值，A^*是A的共轭转置
$$

有人觉得这个定义不太好，因为这个定义是针对特征值的，而不是针对矩阵的；但从理解角度而言这个定义可能确实更好，


```python
# numpy.linalg.norm()函数用于计算矩阵或向量范数, 默认为2范数, 也可以指定ord参数为1范数或无穷范数
# ord = 1, 1范数, 各列绝对值之和的最大值
# ord = np.inf, 无穷范数, 各行绝对值之和的最大值
# ord = -np.inf, 负无穷范数, 各行绝对值之和的最小值
# ord = None, 矩阵范数, 矩阵的最大奇异值
# ord = 'fro', 矩阵范数, 矩阵的Frobenius范数
# ord = 'nuc', 矩阵范数, 矩阵的谱范数

# 测试范数
import numpy as np

# 向量
a = np.array([1, 2, 3])
# 将a转换为矩阵
a1 = np.array([[1, 2, 3]])
# 矩阵
A = np.array([
    [1, 2, 3], 
    [4, 5, 6]
    ])

# 1范数
print(np.linalg.norm(a, ord=1)) # 6.0 如果a多加个中括号会怎么样？
print(np.linalg.norm(a1, ord=1)) # 3.0 这个完全可以说明了其实只有一个中括号时，numpy默认这就是一个列向量
print(np.linalg.norm(A, ord=1, axis=0)) # [5. 7. 9.]

# 无穷范数
print(np.linalg.norm(a, ord=np.inf)) # 3.0 
print(np.linalg.norm(a1, ord=np.inf)) # 6.0
print(np.linalg.norm(A, ord=np.inf)) # 15.0

# 负无穷范数
print(np.linalg.norm(a, ord=-np.inf)) # 1.0
print(np.linalg.norm(a1, ord=-np.inf)) # 6.0
print(np.linalg.norm(A, ord=-np.inf)) # 6.0

# Frobenius范数
print(np.linalg.norm(A, ord='fro')) # 9.539392014169456

# 谱范数
print(np.linalg.norm(A, ord='nuc')) # 10.280901636369205

# 2范数
print(np.linalg.norm(a)) # 3.7416573867739413
print(np.linalg.norm(a1)) # 3.7416573867739413
print(np.linalg.norm(A)) # 9.539392014169456 矩阵范数默认为Frobenius范数

# axis参数
# axis = 0, 按列计算
# axis = 1, 按行计算
# axis = None, 计算整个矩阵范数
print(np.linalg.norm(A, ord=1, axis=0)) # [5. 7. 9.] 按列计算1范数
print(np.linalg.norm(A, ord=1, axis=1)) # [ 6. 15.] 按行计算1范数
print(np.linalg.norm(A, axis=0)) # [4.12310563 5.38516481 6.70820393] 按列计算Frobenius范数
print(np.linalg.norm(A, axis=1)) # [3.74165739 8.77496439] 按行计算Frobenius范数
```

    6.0
    3.0
    [5. 7. 9.]
    3.0
    6.0
    15.0
    1.0
    6.0
    6.0
    9.539392014169456
    10.280901636369205
    3.7416573867739413
    3.7416573867739413
    9.539392014169456
    [5. 7. 9.]
    [ 6. 15.]
    [4.12310563 5.38516481 6.70820393]
    [3.74165739 8.77496439]
    

## 赛德尔迭代法

### 数学推理

赛德尔迭代法是一种迭代法，又叫高斯赛德尔迭代法。

对于方程组$AX=B$，我们可以将其中的A分解为$A = L + D + U$，其中L是下三角矩阵，D是对角矩阵，U是上三角矩阵, 如下所示

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \cdots & a_{nn}
\end{bmatrix} =
\begin{bmatrix}
0 & 0 & \cdots & 0 \\
l_{21} & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
l_{n1} & l_{n2} & \cdots & 0
\end{bmatrix} 
+
\begin{bmatrix}
d_{11} & 0 & \cdots & 0 \\
0 & d_{22} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & d_{nn}
\end{bmatrix}
+
\begin{bmatrix}
0 & u_{12} & \cdots & u_{1n} \\
0 & 0 & \cdots & u_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 0
\end{bmatrix}
$$

其中$l_{ij} = a_{ij}, i>j$，$d_{ii} = a_{ii}, i=1,2,\cdots,n$，$u_{ij} = a_{ij}, i<j$。

则我们可以将方程组$AX=B$转化为$X^{(k+1)} = -D^{-1}LX^{(k+1)} - D^{-1}UX^{(k)} + D^{-1}B$，其中$X^{(k)}$表示第k次迭代的结果。


### 代码实现


```python
import numpy as np

def seidel(A, b, x0, eps=1e-6, max_iter=10):
    """
    高斯赛德尔迭代法
    :param A: 系数矩阵
    :param b: 常数向量
    :param x0: 初始值
    :param tol: 容许误差
    :param max_iter: 最大迭代次数
    :return: 迭代结果
    """
    n = len(b)
    x = x0
    for k in range(max_iter):
        print('第{}次迭代'.format(k))
        print('X = {}'.format(x))
        x_new = np.zeros(n)
        for i in range(n):
            x_new[i] = (b[i] - np.dot(A[i, :i], x_new[:i]) - np.dot(A[i, i+1:], x[i+1:])) / A[i, i]
            print('x = {}'.format(x_new))
        if np.linalg.norm(x_new - x) < eps:
            print('迭代次数：', k)
            break
        x = x_new
    return x

A = np.array([
    [20, 2, 3], 
    [1, 8, 1], 
    [2, -3, 15]
    ])
B = np.array([
    24,
    12, 
    30
    ])
x0 = np.array([0, 0, 0])
x = seidel(A, B, x0)
print(x)
    
```

    第0次迭代
    X = [0 0 0]
    x = [1.2 0.  0. ]
    x = [1.2  1.35 0.  ]
    x = [1.2  1.35 2.11]
    第1次迭代
    X = [1.2  1.35 2.11]
    x = [0.7485 0.     0.    ]
    x = [0.7485    1.1426875 0.       ]
    x = [0.7485    1.1426875 2.1287375]
    第2次迭代
    X = [0.7485    1.1426875 2.1287375]
    x = [0.76642062 0.         0.        ]
    x = [0.76642062 1.13810523 0.        ]
    x = [0.76642062 1.13810523 2.12543163]
    第3次迭代
    X = [0.76642062 1.13810523 2.12543163]
    x = [0.76737473 0.         0.        ]
    x = [0.76737473 1.1383992  0.        ]
    x = [0.76737473 1.1383992  2.12536321]
    第4次迭代
    X = [0.76737473 1.1383992  2.12536321]
    x = [0.7673556 0.        0.       ]
    x = [0.7673556  1.13841015 0.        ]
    x = [0.7673556  1.13841015 2.12536795]
    第5次迭代
    X = [0.7673556  1.13841015 2.12536795]
    x = [0.76735379 0.         0.        ]
    x = [0.76735379 1.13840978 0.        ]
    x = [0.76735379 1.13840978 2.12536812]
    第6次迭代
    X = [0.76735379 1.13840978 2.12536812]
    x = [0.7673538 0.        0.       ]
    x = [0.7673538  1.13840976 0.        ]
    x = [0.7673538  1.13840976 2.12536811]
    迭代次数： 6
    [0.76735379 1.13840978 2.12536812]
    


```python
# 对比一下雅可比迭代法
jacobi(A, B, x0) # 迭代次数8次
```

    第0次迭代
    X = [0 0 0]
    第1次迭代
    X = [1.2 1.5 2. ]
    第2次迭代
    X = [0.75 1.1  2.14]
    第3次迭代
    X = [0.769   1.13875 2.12   ]
    第4次迭代
    X = [0.768125   1.138875   2.12521667]
    第5次迭代
    X = [0.76733    1.13833229 2.12535833]
    第6次迭代
    X = [0.76736302 1.13841396 2.12535579]
    第7次迭代
    X = [0.76735524 1.13841015 2.12536772]
    第8次迭代
    X = [0.76735383 1.13840963 2.125368  ]
    迭代次数： 8
    




    array([0.76735384, 1.13840977, 2.12536808])



## 松弛迭代法


### 数学推理

松弛迭代法就是在残差的基础上加入了松弛因子$\omega$，

$$
X^{(k+1)} = X^{(k)} + \omega R^{(k)}
$$


### 简单迭代法下的逐次松弛法

对于雅可比迭代法的迭代公式, 我们可以这样写

$$
\begin{aligned}
X^{(k+1)} &= (I -  D^{-1}A)X^{(k)} +  D^{-1}B \\
&= X^{(k)} +  D^{-1}(B - AX^{(k)})
\end{aligned}
$$

我们把$D^{-1}(B - AX^{(k)})$称为残差，因此我们可以将迭代公式写成

$$
X^{(k+1)} = X^{(k)} + R^{(k)}
$$

引入松弛因子, 我们可以得到

$$
\begin{aligned}
X^{(k+1)} &= X^{(k)} + \omega R^{(k)} \\
&= X^{(k)} + \omega D^{-1}(B - AX^{(k)})
\end{aligned}
$$

#### 代码实现



```python

def jacobi_relax(A, b, x0, omega, max_iter=100, tol=1e-6):
    n = len(b)
    x = x0
    for i in range(max_iter):
        x_new = x + omega * np.linalg.inv(np.diag(np.diag(A))).dot(b - A.dot(x))
        if np.linalg.norm(x_new - x) < tol:
            print('迭代次数：', i)
            break
        x = x_new
    return x

```


```python
jacobi_relax(A, B, x0, 0.5) # 迭代次数24次
```

    迭代次数： 24
    




    array([0.76735483, 1.13840975, 2.12536716])




```python
jacobi_relax(A, B, x0, 1) # 迭代次数8次
```

    迭代次数： 8
    




    array([0.76735383, 1.13840963, 2.125368  ])




```python
jacobi_relax(A, B, x0, 1.25) # 迭代次数16次
```

    迭代次数： 16
    




    array([0.76735393, 1.13840993, 2.12536793])



### 高斯赛德尔迭代法下的逐次松弛法

对于高斯赛德尔迭代法的迭代公式, 我们可以这样写

$$
\begin{aligned}
X^{(k+1)} &= -D^{-1}(L)X^{(k+1)} - D^{-1}(U)X^{(k)} + D^{-1}B \\
&= D^{-1}(B - LX^{(k+1)} - UX^{(k)})
\end{aligned}
$$

对于单个元素，我们可以这样写

$$
\begin{aligned}
x_{i}^{(k+1)} &= d_{ii}^{-1}(b_{i} - \sum_{j=1}^{i-1}l_{ij}x_{j}^{(k+1)} - \sum_{j=i+1}^{n}u_{ij}x_{j}^{(k)}) \\
&= \frac{1}{a_{ii}}(b_{i} - \sum_{j=1}^{i-1}l_{ij}x_{j}^{(k+1)} - \sum_{j=i+1}^{n}u_{ij}x_{j}^{(k)})
\\ \Rightarrow
x_{i}^{(k+1)} &= x_{i}^{(k)} + \frac{1}{a_{ii}}(b_{i} - \sum_{j=1}^{i-1}l_{ij}x_{j}^{(k+1)} -a_{ii}x_{i}^{(k)}- \sum_{j=i+1}^{n}u_{ij}x_{j}^{(k)})
\\ 引入松弛因子\omega \Rightarrow
x_{i}^{(k+1)} &= x_{i}^{(k)} + \frac{\omega}{a_{ii}}(b_{i} - \sum_{j=1}^{i-1}l_{ij}x_{j}^{(k+1)} -a_{ii}x_{i}^{(k)}- \sum_{j=i+1}^{n}u_{ij}x_{j}^{(k)})
\\ \Rightarrow
x_{i}^{(k+1)} &= x_{i}^{(k)} + \omega (\hat{x}_i - x_i^{(k)})  = (1-\omega)x_i^{(k)} + \omega \hat{x}_i
\\ \text{其中} \hat{x}_i &= \frac{1}{a_{ii}}(b_{i} - \sum_{j=1}^{i-1}l_{ij}x_{j}^{(k+1)} - \sum_{j=i+1}^{n}u_{ij}x_{j}^{(k)}) 
\end{aligned}
$$



#### 代码实现



```python
def gauss_seidel_relaxation(A, b, x0, omega, max_iter=100, tol=1e-6):
    """
    高斯赛德尔迭代法下的逐次松弛法
    :param A: 系数矩阵
    :param b: 常数向量
    :param x0: 初始值
    :param omega: 松弛因子
    :param max_iter: 最大迭代次数
    :param tol: 精度
    :return: 迭代结果
    """
    n = len(b)
    x = x0
    for i in range(max_iter):
        x_new = np.zeros(n)
        for j in range(n):
            x_new[j] = (1 - omega) * x[j] + omega * (b[j] - np.dot(A[j, :j], x_new[:j]) - np.dot(A[j, j + 1:], x[j + 1:])) / A[j, j]
        if np.linalg.norm(x_new - x) < tol:
            print('迭代次数：', i)
            break
        x = x_new
    return x
```


```python
gauss_seidel_relaxation(A, B, x0, 0.5) # 迭代次数22次
```

    迭代次数： 22
    




    array([0.76735551, 1.13841036, 2.125367  ])




```python
gauss_seidel_relaxation(A, B, x0, 1) # 迭代次数6次
```

    迭代次数： 6
    




    array([0.76735379, 1.13840978, 2.12536812])




```python
gauss_seidel_relaxation(A, B, x0, 1.25) # 迭代次数13次
```

    迭代次数： 13
    




    array([0.76735303, 1.13840981, 2.12536829])



> 注: 对omega的选择, 有一些经验值, 例如当$\omega > 2$时, 逐次松弛法会发散, 因此, 一般选择$\omega \in (0, 2)$ , 当$\omega = 1$时, 逐次松弛法退化为简单迭代法或者赛德尔;
> 数值分析课本证明了$\omega \in (0, 2)$是必要条件, 但是没有证明充分条件;

## 例2

<!-- 验证课本例题 -->

用赛德尔迭代法下的松弛迭代解下列方程组

$$
\begin{cases}
4x_1 + 3x_3 = 24 \\
3x_1 + 4x_2 - x_3 = 30 \\
-x_2 + 4x_3 = -24
\end{cases}
$$


### 代码实现


```python
import numpy as np

def gauss_seidel(A, b, x0, omega=1, tol=1e-7, max_iter=100):
    n = A.shape[0]
    x = x0.copy()
    for i in range(max_iter):
        print('X = {}'.format(x))
        for j in range(n):
            x[j] = (1 - omega) * x[j] + omega * (b[j] - A[j, :j] @ x[:j] - A[j, j+1:] @ x[j+1:]) / A[j, j]
            # x_new[j] = (1 - omega) * x[j] + omega * (b[j] - np.dot(A[j, :j], x_new[:j]) - np.dot(A[j, j + 1:], x[j + 1:])) / A[j, j] 
        if np.linalg.norm(A @ x - b) < tol:
            print('迭代次数：', i)
            break
    return x

```


```python
A = np.array([[4, 0, 3], [3, 4, -1], [0, -1, 4]])
b = np.array([24, 30, -24])
x0 = np.array([1,1,1], dtype=float)
x = gauss_seidel(A, b, x0, omega=1.25)
print("迭代结果：", x)
print("误差：", np.linalg.norm(A @ x - b))
```

    X = [1. 1. 1.]
    X = [ 6.3125      3.51953125 -6.65014648]
    X = [12.15638733 -4.97966671 -7.39360923]
    X = [11.39241182 -2.37097228 -6.39252653]
    X = [10.64489067 -2.00950647 -6.52983914]
    X = [10.96050153 -2.43866829 -6.62962406]
    X = [10.97514717 -2.37629092 -6.5851849 ]
    X = [10.92982405 -2.3355076  -6.5835499 ]
    X = [10.93962202 -2.35437809 -6.58985568]
    X = [10.94308419 -2.35487681 -6.58843508]
    X = [10.94088684 -2.35224818 -6.58796878]
    X = [10.94099902 -2.35286479 -6.58827805]
    X = [10.94126092 -2.3530528  -6.58825949]
    X = [10.94117804 -2.3529223  -6.58822335]
    X = [10.94116488 -2.35293129 -6.58823519]
    X = [10.94117927 -2.35294624 -6.5882369 ]
    X = [10.94117728 -2.35294117 -6.58823489]
    X = [10.94117589 -2.35294051 -6.58823519]
    X = [10.94117651 -2.35294135 -6.58823538]
    X = [10.94117654 -2.35294122 -6.58823529]
    迭代次数： 19
    迭代结果： [10.94117645 -2.35294114 -6.58823529]
    误差： 8.793428993778205e-08
    


```python
np.linalg.solve(A, b)
```




    array([10.94117647, -2.35294118, -6.58823529])




```python
gauss_seidel(A, b, x0, omega=1)
```

    X = [1. 1. 1.]
    X = [ 5.25      3.8125   -5.046875]
    X = [ 9.78515625 -1.10058594 -6.27514648]
    X = [10.70635986 -2.09855652 -6.52463913]
    X = [10.89347935 -2.30126929 -6.57531732]
    X = [10.93148799 -2.34244533 -6.58561133]
    X = [10.9392085  -2.35080921 -6.5877023 ]
    X = [10.94077673 -2.35250812 -6.58812703]
    X = [10.94109527 -2.35285321 -6.5882133 ]
    X = [10.94115998 -2.35292331 -6.58823083]
    X = [10.94117312 -2.35293755 -6.58823439]
    X = [10.94117579 -2.35294044 -6.58823511]
    X = [10.94117633 -2.35294103 -6.58823526]
    迭代次数： 12
    




    array([10.94117644, -2.35294115, -6.58823529])



### 小结

个人感觉课本上的习题的解释应该是错了，就比如例2，课本上的解释是

<!-- 引用 -->

> 若要迭代结果精确到7位小数，高斯赛德尔迭代法（选取$\omega = 1$）需要34次迭代运算，而松弛迭代法（选取$\omega = 1.25$）只需要迭代14次，

而实际上，上面的代码的结果也看到了，实际上是松弛迭代法需要19次迭代，而高斯赛德尔迭代法只需要12次迭代，松弛迭代法的收敛速度也不见得比高斯赛德尔迭代法快。

> 个人感觉而已，也有可能是代码有问题，欢迎指正。
