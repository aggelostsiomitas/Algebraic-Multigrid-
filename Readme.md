This codes solves the Ax=b using the AMG Multigrid method 
The code contains the below set of functions:

$\bullet$
**set_intersection function:**
Ths function compute the intersection between two vectors $A\cap B$ \
$\textbf{Inputs}$: \
$\hspace{5cm} \bullet subspace A $ \
$\hspace{5cm} \bullet subspace B $ 

$\bullet$
**set_difference function**
This function computer the difference between two subspaces $A\setminus B$ \
$\textbf{Inputs}$: \
$\hspace{5cm} \bullet subspace A $ \
$\hspace{5cm} \bullet subspace B $ 

$\bullet$
**find function**
This function searches the rows of the matrix and finds the non zero element's columns and values 
It returns the columns and values as a pair \
$\textbf{Inputs}$: \
            $\hspace{5cm} \bullet$ coefficient matrix A \
            $\hspace{5cm}\bullet$ row

$\bullet$
**connecton function**
This function computer the connection between the nodes of the matrix A
The function has as input the coefficient matrix  A and returns the S,ST coefficient matrices 

$\textbf{Inputs}$: \
            $\hspace{5cm} \bullet $ coefficient matrix A 

$\textbf{Outputs}$: \
            $\hspace{5cm} \bullet $ coefficient matrix S that will hold the strong variables that influence node i\
            $\hspace{5cm} \bullet $ coefficient matrix ST that will hold the variables that strongly get influenced by node i

The mathematical formula for the coefficient matrix S is : $S_{i}=-a_{ij}>=\theta \cdot max(-a_{ik})$.
The coefficient matrices S,ST will be computer by searching the matrix A and finding the non zero elements's cols,values using the find function.A vector strength will hold the negative values of the values found computing the term $-a_{ij}$ and a variable alpha will hold the maximum value computing the term $-a_{ik}$.The theta variable exists in the variables of the class AMGSOLVER  The columns that agree with the formula : $-a_{ij}>=\theta \cdot max(-a_{ik})$  will be the values of S and for every $S_i$ ST will hold the i values: $ST_j=i \quad \forall j \in  S_i$

$\bullet$
**coarsening function**
This function will compute the coarsening between the coarse levels and fine levels.
$\textbf{Inputs}$: \
            $\hspace{5cm} \bullet$ coefficient matrix S \
            $\hspace{5cm}\bullet$ coefficient matrix ST

$\textbf{Output}$: \
            $\hspace{5cm} \bullet$ coefficient vector C that willl hold the nodes that willbecome coarse  \
            $\hspace{5cm}\bullet$ coefficient vector F that will hold the nodes that will become fine 

The coarsening will happen with the use of the Ruge-Stuben algorithm that states that: 

$$\begin{aligned} 
& \textbf{Coarsening algorithm} \\
& \textbf{Inputs: } \text{coefficient matrix }S ,\text{coefficient matrix }ST \\ 
& \textbf{Output: }\text{coefficient vector }C , \text{coefficient vector }F \\
&1: \lambda_i=\Vert{} ST_i\Vert{} \\
&2: \text{Keep track of the state of the nodes : 0=undecided,1=coarse,-1=fine in a state vector}\\
&3: \text{Initialize :} C=\emptyset , F=\emptyset \\
&\text{Start first sweep} \\
&4: \text{while } \text{state } \neq 0 \textbf{ do: } \\ 
&5:\quad \quad \textbf{for } i \in \text{state} \textbf{ do: }\\
&\quad \quad \quad \quad   i=\arg \max_k(\lambda_k) \\
&\quad \quad \quad  \quad \text{Set } \text{state(i)}=1 \\
&6:\quad \quad \quad \forall j\in ST_i \textbf{ If } \text{state(j)=0} \textbf{ do: } \\
&\quad \quad \quad \quad \quad \quad \text{state(j)}=-1 \\
&\quad \quad \quad \quad \quad \quad \forall k \in ST_j \textbf{ If } \text{state(k)=0} \textbf{ do: }\\
&7: \quad \quad \quad \quad \quad \quad \quad \lambda_k=\lambda_k +1\\
&\quad \quad \quad \quad \quad \quad \textbf{ End do} \\
&\quad \quad \quad \quad  \textbf{ End do} \\
&\quad \quad \quad   \textbf{ End do} \\
&\quad \quad    \textbf{ End do} \\
& \\
&\text{Start second sweep }\\
&8: \forall i \in \{\text{state} \mid \text{state(i)=1}\} \textbf{ do: }\\
&9:\quad \forall j \in S_i \cap \{state=-1\} \textbf{ do: }\\
& \quad \quad  \quad \textbf{ If } S_i \cap \{S_j \cap \{state=1\} \} =\emptyset \textbf{ do: }\\
& \quad \quad  \quad \quad \text{Set state(i)=1}\\
& \quad \quad  \quad \textbf{End do} \\
& \quad \quad  \textbf{End do} \\
& \quad \textbf{End do} \\
&10: C=C\cap \{state=1\} , F=F\cap \{state=-1\}
\end{aligned}$$

$\bullet$
**interpolate function**
The interpolation function will compute the interpolation between the coarse and fine points \
$\text{Inputs: }$ \
            $\hspace{5cm} \bullet$ coefficient matrix A \
            $\hspace{5cm}\bullet$ coefficient matrix S \
            $\hspace{5cm}\bullet$ coefficient vector C \
            $\hspace{5cm}\bullet$ coefficient vector F 

$\text{Outputs: }$  
            $\hspace{5cm} \bullet$ coefficient matrix P 

The interpolation function will compute the interpolation by computing the necessary weights for the nodes. The matreix P will have dimension nxnc where n is the size of the matrix A and nc the size of the vector C that was computed in the coarsening function 
The function will initialize the map vector with dimensions nx1 that will map the global indexes to P indexes 
$$\forall k \in nc : map(C(k))=k$$
For the coarse points: 
$$P(C(k),k)=1  \quad \forall k \in nc$$

For the computation of the weights the following algorithm will be used 
for the computation of \
 $P_(i,k)=\frac{A(i,j)+\sum_{m\in F_i^s }\left(\frac{A(i,m)\cdot A(m,j)}{\sum_{k \in C_i}A(m,k)} \right)}{A(i,i)+\sum_{m\in F_i^w} A(i,m)} \quad \forall j \in C_i $ 

The algorithm is shown below: 

$$\begin{aligned} 
&1  \forall i \in F \textbf{ do: } \\
&2 \quad \text{Compute: } F_i =S_i \cap F ,\quad  C_i =S_i \cap C \\
&3 \quad \text{Compute the weak nodes: } F_{i}^w=F-F_i^s \\
&4 \quad \text{Compute :} down=A(i,i)+\sum_{m \in F_i^w} A_{im} \\
&5 \quad \forall j\in Ci \textbf{ do: } \\
&6 \quad  \quad \text{Compute the term: } up=A(i,j) \\
&7 \quad  \quad \forall m\in F_i^s \textbf{ do: } \\
&8 \quad  \quad \quad up=up+\frac{A(i,m)*A(m,j)}{\sum_{k \in C_i}A(m,k)} \\
& \quad \quad \quad \textbf{End do} \\
&9 \quad \quad P(i,map(j))=-up.down \\
& \quad \quad \textbf{End do} \\
& \quad  \textbf{End do} \\
\end{aligned}$$

$\bullet$
**Galrking function**
$\text{Inputs: }$ \
            $\hspace{5cm} \bullet$ coefficient matrix A \
            $\hspace{5cm}\bullet$ coefficient matrix P\

$\text{Outputs: }$  
            $\hspace{5cm} \bullet$ coefficient matrix $A_c$ \
The galerkin function will compute the mathematical computation of $$A_c=P^T\cdot A\cdot P$$

