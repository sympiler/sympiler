//
// Created by kazem on 6/14/17.
//
#include "NumericalUtils.h"

/* depth-first-search of the graph of a matrix, starting at node j */
int dfs (int j, int *Gp, int *Gi, int top, int *xi, int *pstack, const int *pinv)
{
    int i, p, p2, done, jnew, head = 0 ;
    //if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
    //Gp = G->p ; Gi = G->i ;
    xi [0] = j ;                /* initialize the recursion stack */

    while (head >= 0)
    {
        j = xi [head] ;         /* get j from the top of the recursion stack */
        jnew = pinv ? (pinv [j]) : j ;
        if (!CS_MARKED (Gp, j))
        {
            CS_MARK (Gp, j) ;       /* mark node j as visited */
            pstack [head] = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew]) ;
        }
        done = 1 ;                  /* node j done if no unvisited neighbors */
        p2 = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew+1]) ;
        for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
        {
            i = Gi [p] ;            /* consider neighbor node i */
            if (CS_MARKED (Gp, i)) continue ;   /* skip visited node i */
            pstack [head] = p ;     /* pause depth-first search of node j */
            xi [++head] = i ;       /* start dfs at node i */
            done = 0 ;              /* node j is not done */
            break ;                 /* break, to start dfs (i) */
        }
        if (done)               /* depth-first search at node j is done */
        {
            head-- ;            /* remove j from the recursion stack */
            xi [--top] = j ;    /* and place in the output stack */
        }
    }
    return (top) ;
}

/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
int reach (int n, int* Gp, int *Gi, int* Bp, int* Bi, int k, int *pruneVal, const int *pinv)
{
    int p, top, setSize=0 ;
    int *xi = new int[2*n]();
    if (!Gp || !Gi || !Bp || !Bi || !xi) return (-1) ;    /* check inputs */
    //n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
    top = n ;
    for (p = Bp [k] ; p < Bp [k+1] ; p++)
    {
        if (!CS_MARKED (Gp, Bi [p]))    /* start a dfs at unmarked node i */
        {
            top = dfs (Bi [p], Gp, Gi, top, xi, xi+n, pinv) ;
        }
    }
    for (p = top ; p < n ; p++) CS_MARK (Gp, xi [p]) ;  /* restore G */
    for (int i = top; i < n; ++i) {
        pruneVal[setSize++];
    }
    delete []xi;
    return (--setSize) ;
}

int reach_sn (int n, int* Gp, int *Gi, int* Bp, int* Bi,
              int k, int *PBset, const int *pinv, int sn,
              int *col2sup ){
    int p, top, PBsize=0 ;
    double analysis;
    bool *checked = new bool[sn]();
    int *xi = new int[2*n]();
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    start = std::chrono::system_clock::now();
    if (!Gp || !Gi || !Bp || !Bi || !xi) return (-1) ;    /* check inputs */
    //n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
    top = n ;
    for (p = Bp [k] ; p < Bp [k+1] ; p++)
    {
        if (!CS_MARKED (Gp, Bi [p]))    /* start a dfs at unmarked node i */
        {
            top = dfs (Bi [p], Gp, Gi, top, xi, xi+n, pinv) ;
        }
    }
    for (p = top ; p < n ; p++) CS_MARK (Gp, xi [p]) ;  /* restore G */
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    analysis=elapsed_seconds.count();
    for (int i = top; i < n; ++i) {
        if(!checked[col2sup[xi[i]]]){
            checked[col2sup[xi[i]]] = true;
            PBset[PBsize++] = col2sup[xi[i]];
        }
    }
    delete []checked;
    delete []xi;
    return (PBsize) ;
}

int* superNodeDetection(int *Parent, int *ColCount, int n,
                        int *col2sup, int &bsSize, int &avgBSize){
    int *Wi = new int[n]();
    int *blockSet = new int[n]();
    int blockSetSize=0, j=0, parent;
    avgBSize=0;
    /* ---------------------------------------------------------------------- */
    /* find the fundamental supernodes */
    /* ---------------------------------------------------------------------- */

    /* count the number of children of each node, using Wi [ */
    for (j = 0 ; j < n ; j++)
    {
        parent = Parent [j] ;
        if (parent != EMPTY)
        {
            Wi [parent]++ ;
        }
    }

    //blockSet = Head ;  /* use Head [0..bsSize] as workspace for blockSet list ( */

    /* column 0 always starts a new supernode */
    bsSize = (n == 0) ? 0 : 1 ;	/* number of fundamental supernodes */
    blockSet [0] = 0 ;

    for (j = 1 ; j < n ; j++)
    {
        /* check if j starts new supernode, or in the same supernode as j-1 */
        if (Parent [j-1] != j	    /* parent of j-1 is not j */
            || (ColCount [j-1] != ColCount [j] + 1) /* j-1 not subset of j*/
            || Wi [j] > 1	    /* j has more than one child */
                )
        {
            /* j is the leading node of a supernode */
            blockSet [bsSize++] = j ;
        }
    }
    blockSet [bsSize] = n ;

    /* contents of Wi no longer needed for child count ] */

    /* ---------------------------------------------------------------------- */
    /* find the mapping of fundamental nodes to supernodes */
    /* ---------------------------------------------------------------------- */

    //col2sup = Wj ;	/* use Wj as workspace for col2sup [ */
    //int *col2sup = new int[bsSize]();
    /* col2sup [k] = s if column k is contained in supernode s */
    for (int s = 0 ; s < bsSize ; s++)
    {
        for (int k = blockSet [s] ; k < blockSet [s+1] ; k++)
        {
            col2sup [k] = s ;
        }
        avgBSize+=(blockSet [s+1] - blockSet [s]);
    }
    avgBSize/=bsSize;
    return blockSet;
}

int *etree (int n, int* Ap, int* Ai, int ata){
    /* compute the etree of A (using triu(A),
     * or A'A without forming A'A
     * n: the matrix size
     * Ap: column pointer
     * Ai: row index
     * ata: A'A or A
     * */
    int i, k, p, m, inext, *w, *parent, *ancestor, *prev ;
    if (n<0 || Ap == NULL || Ai == NULL) //check inputs
        return 0;
    m = n;
    parent = new int[n]; //result allocation
    w = new int[n + (ata ? m : 0)]; // get workspace
    if (w == NULL || parent == NULL)
        return 0;

    ancestor = w ; prev = w + n ;
    if (ata) for (i = 0 ; i < m ; i++) prev [i] = -1 ;
    for (k = 0 ; k < n ; k++)
    {
        parent [k] = -1 ;                   /* node k has no parent yet */
        ancestor [k] = -1 ;                 /* nor does k have an ancestor */
        for (p = Ap [k] ; p < Ap [k+1] ; p++)
        {
            i = ata ? (prev [Ai [p]]) : (Ai [p]) ;
            for ( ; i != -1 && i < k ; i = inext)   /* traverse from i to k */
            {
                inext = ancestor [i] ;              /* inext = ancestor of i */
                ancestor [i] = k ;                  /* path compression */
                if (inext == -1) parent [i] = k ;   /* no anc., parent is k */
            }
            if (ata) prev [Ai [p]] = k ;
        }
    }
    delete []w;
    return parent;
}

/* depth-first search and postorder of a tree rooted at node j */
int tdfs(int j, int k, int *head, const int *next, int *post, int *stack)
{
    int i, p, top = 0 ;
    if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
    stack [0] = j ;                 /* place j on the stack */
    while (top >= 0)                /* while (stack is not empty) */
    {
        p = stack [top] ;           /* p = top of stack */
        i = head [p] ;              /* i = youngest child of p */
        if (i == -1)
        {
            top-- ;                 /* p has no unordered children left */
            post [k++] = p ;        /* node p is the kth postordered node */
        }
        else
        {
            head [p] = next [i] ;   /* remove i from children of p */
            stack [++top] = i ;     /* start dfs on child node i */
        }
    }
    return (k) ;
}

int *postOrder (const int *parent, int n)
{
/* post order a forest
 * Obtained from CSparse library
 * */
    int j, k = 0, *post, *w, *head, *next, *stack ;
    if (!parent) return (NULL) ;                        /* check inputs */
    //post = cs_malloc (n, sizeof (csi)) ;                /* allocate result */
    //w = cs_malloc (3*n, sizeof (csi)) ;                 /* get workspace */
    post = new int[n];
    w = new int[3*n];
    //if (!w || !post) return (cs_idone (post, NULL, w, 0)) ;
    if (!w || !post)
        return NULL;
    head = w ; next = w + n ; stack = w + 2*n ;
    for (j = 0 ; j < n ; j++) head [j] = -1 ;           /* empty linked lists */
    for (j = n-1 ; j >= 0 ; j--)            /* traverse nodes in reverse order*/
    {
        if (parent [j] == -1) continue ;    /* j is a root */
        next [j] = head [parent [j]] ;      /* add j to list of its parent */
        head [parent [j]] = j ;
    }
    for (j = 0 ; j < n ; j++)
    {
        if (parent [j] != -1) continue ;    /* skip j if it is not a root */
        k = tdfs(j, k, head, next, post, stack) ;
    }
    // return (cs_idone (post, NULL, w, 1)) ;  /* success; free w, return post */
    delete []w;
    return post;
}


double cumsum (int *p, int *c, int n)
{
    int i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid csi overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

int transpose (int n, int m, int *Ap, int *Ai, double *Ax, int values,
               int *Cp, int *Ci, double *Cx)
{
    int p, q, j, *w ;
    //double *Cx;
    //cs *C ;
    if (!Ai || !Ap) return ((int)NULL) ;    /* check inputs */
    //m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    //C = cs_spalloc (n, m, Ap [n], values && Ax, 0) ;       /* allocate result */
    //Cp = new int[n+1]; Ci = new int[Ap[n]];
    //if(values) Cx = new double[Ap[n]];
    //else Cx=NULL;
    //w = cs_calloc (m, sizeof (csi)) ;                      /* get workspace */
    w = new int[m]();
    if (!Cp || !w) return -1 ;       /* out of memory */
    //Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    delete []w;
    return 1 ;  /* success;  */
}

int transposePattern (int n, int m, int *Ap, int *Ai, int values,
               int *Cp, int *Ci, double *Cx)
{
    int p, q, j, *w ;
    //double *Cx;
    //cs *C ;
    if (!Ai || !Ap) return ((int)NULL) ;    /* check inputs */
    //m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    //C = cs_spalloc (n, m, Ap [n], values && Ax, 0) ;       /* allocate result */
    //Cp = new int[n+1]; Ci = new int[Ap[n]];
    //if(values) Cx = new double[Ap[n]];
    //else Cx=NULL;
    //w = cs_calloc (m, sizeof (csi)) ;                      /* get workspace */
    w = new int[m]();
    if (!Cp || !w) return -1 ;       /* out of memory */
    //Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            //if (Cx) Cx [q] = Ax [p] ;
        }
    }
    delete []w;
    return 1 ;  /* success;  */
}

/* consider A(i,j), node j in ith row subtree and return lca(jprev,j) */
int leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
          int *ancestor, int *jleaf)
{
    int q, s, sparent, jprev ;
    if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1) ;
    *jleaf = 0 ;
    if (i <= j || first [j] <= maxfirst [i]) return (-1) ;  /* j not a leaf */
    maxfirst [i] = first [j] ;      /* update max first[j] seen so far */
    jprev = prevleaf [i] ;          /* jprev = previous leaf of ith subtree */
    prevleaf [i] = j ;
    *jleaf = (jprev == -1) ? 1: 2 ; /* j is first or subsequent leaf */
    if (*jleaf == 1) return (i) ;   /* if 1st leaf, q = root of ith subtree */
    for (q = jprev ; q != ancestor [q] ; q = ancestor [q]) ;
    for (s = jprev ; s != q ; s = sparent)
    {
        sparent = ancestor [s] ;    /* path compression */
        ancestor [s] = q ;
    }
    return (q) ;                    /* q = least common ancester (jprev,j) */
}


void init_ata (int n, int m, int *ATp, int *ATi, const int *post, int *w, int **head, int **next)
{
    int i, k, p ;
    *head = w+4*n, *next = w+5*n+1 ;
    for (k = 0 ; k < n ; k++) w [post [k]] = k ;    /* invert post */
    for (i = 0 ; i < m ; i++)
    {
        for (k = n, p = ATp[i] ; p < ATp[i+1] ; p++) k = CS_MIN (k, w [ATi[p]]);
        (*next) [i] = (*head) [k] ;     /* place row i in linked list k */
        (*head) [k] = i ;
    }
}
int *counts (int n, int m, int *Ap, int *Ai, const int *parent, const int *post, int ata)
{
    int i, j, k, J, s, p, q, jleaf, *ATp, *ATi, *maxfirst, *prevleaf,
            *ancestor, *head = NULL, *next = NULL, *colcount, *w, *first, *delta ;
    //cs *AT ;
    if (!Ap || !parent || !post) return (NULL) ;    /* check inputs */
    //m = A->m ; n = A->n ;
    s = 4*n + (ata ? (n+m+1) : 0) ;
    //delta = colcount = cs_malloc (n, sizeof (int)) ;    /* allocate result */
    delta = colcount = new int[n];
    w = new int[s];                   /* get workspace */
    ATp = new int[n+1]; ATi = new int[Ap[n]];
    transposePattern (n,m,Ap,Ai,0,ATp,ATi,NULL);                          /* AT = A' */
    if (!ATp || !ATi || !colcount || !w) return (colcount) ;
    ancestor = w ; maxfirst = w+n ; prevleaf = w+2*n ; first = w+3*n ;
    for (k = 0 ; k < s ; k++) w [k] = -1 ;      /* clear workspace w [0..s-1] */
    for (k = 0 ; k < n ; k++)                   /* find first [j] */
    {
        j = post [k] ;
        delta [j] = (first [j] == -1) ? 1 : 0 ;  /* delta[j]=1 if j is a leaf */
        for ( ; j != -1 && first [j] == -1 ; j = parent [j]) first [j] = k ;
    }
//    ATp = AT->p ; ATi = AT->i ;
    if (ata) init_ata (n, m, ATp, ATi, post, w, &head, &next) ;
    for (i = 0 ; i < n ; i++) ancestor [i] = i ; /* each node in its own set */
    for (k = 0 ; k < n ; k++)
    {
        j = post [k] ;          /* j is the kth node in postordered etree */
        if (parent [j] != -1) delta [parent [j]]-- ;    /* j is not a root */
        for (J = HEAD (k,j) ; J != -1 ; J = NEXT (J))   /* J=j for LL'=A case */
        {
            for (p = ATp [J] ; p < ATp [J+1] ; p++)
            {
                i = ATi [p] ;
                q = leaf (i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
                if (jleaf >= 1) delta [j]++ ;   /* A(i,j) is in skeleton */
                if (jleaf == 2) delta [q]-- ;   /* account for overlap in q */
            }
        }
        if (parent [j] != -1) ancestor [j] = parent [j] ;
    }
    for (j = 0 ; j < n ; j++)           /* sum up delta's of each child */
    {
        if (parent [j] != -1) colcount [parent [j]] += colcount [j] ;
    }
    delete [] w;
    delete [] ATi;
    delete [] ATp;
    return (colcount) ;    /* success: free workspace */
}


//Reach set for Cholesky

/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
int ereach (int n, int *Ap, int *Ai, int k, const int *parent, int *s, int *w)
{
    int i, p, len, top;
    if (!Ap || !Ai || !parent || !s || !w) return (-1) ;   /* check inputs */
    top = n;
    CS_MARK (w, k) ;                /* mark node k as visited */
    for (p = Ap [k] ; p < Ap [k+1] ; p++)
    {
        i = Ai [p] ;                /* A(i,k) is nonzero */
        if (i > k) continue ;       /* only use upper triangular part of A */
        for (len = 0 ; !CS_MARKED (w,i) ; i = parent [i]) /* traverse up etree*/
        {
            s [len++] = i ;         /* L(k,i) is nonzero */
            CS_MARK (w, i) ;        /* mark i as visited */
        }
        while (len > 0) s [--top] = s [--len] ; /* push path onto stack */
    }
    for (p = top ; p < n ; p++) CS_MARK (w, s [p]) ;    /* unmark all nodes */
    CS_MARK (w, k) ;                /* unmark node k */
    return (top) ;                  /* s [top..n-1] contains pattern of L(k,:)*/
}

/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
int ereach_sn (int n, int *Ap, int *Ai, int col1,int col2, int *col2sup,
               const int *parent, int *s, int *w)
{
    int i, p, len, top;
    if (!Ap || !Ai || !parent || !s || !w) return (-1) ;   /* check inputs */
    top = n;
    for (int k = col1; k < col2; ++k) {
        if(k==col1)
        CS_MARK (w, col2sup[k]) ;                /* mark node k as visited */
        for (p = Ap [k] ; p < Ap [k+1] ; p++){
            i = col2sup[Ai [p]] ;                /* A(i,k) is nonzero block */
            if (Ai [p] > k) continue ;       /* only use upper triangular part of A */
            //if(col2sup[i] == col2sup[Ai[p-1]]) continue; // from the same supernode
            for (len = 0 ; !CS_MARKED (w,i) ; i = parent [i]) /* traverse up etree*/
            {
                s [len++] = i ;         /* L(k,i) is nonzero */
                CS_MARK (w, i) ;        /* mark i as visited */
            }
            while (len > 0) s [--top] = s [--len] ; /* push path onto stack */
        }
    }
    for (p = top ; p < n ; p++) CS_MARK (w, s [p]) ;    /* unmark all nodes */
    //for (int k = col1; k < col2; ++k) {
    CS_MARK (w, col2sup[col1]);                /* unmark node k */
    //}
    return (top) ;                  /* s [top..n-1] contains pattern of L(k,:)*/
}

bool pruning(int n, int* c, int* r, int *ET, int* &prunePtr, int* &pruneSet){
    //Specifying the row patterns of L and pruneSet
    int top;
    int *xi = new int[2*n];
    prunePtr[0]=prunePtr[1]=0;
    for (int colNo = 1; colNo < n; ++colNo) {
        top = ereach(n,c,r,colNo,ET,xi,xi+n);
        prunePtr[colNo+1]=prunePtr[colNo]+(n-top);
        for (int i = top, cnt=0; i < n; ++i, ++cnt) {
            pruneSet[prunePtr[colNo]+cnt]=xi[i];
        }
    }
}