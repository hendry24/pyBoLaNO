from sympy import \
    Add, \
    Mul, \
    Integer, \
    Expr, \
    sympify, \
    gcd
from .operators import \
    is_ladder_contained
from .normal_ordering import \
    normal_order as NO

__all__ = ["expval"]

def get_scalars_and_ops(q):
    """
    q is an input to expval that is not a sum. If the input is a sum, q
    input here must be one of its summands.
    """
    scalars = []
    ops = []

    if not(is_ladder_contained(q)):
        scalars.append(q)
    elif not(isinstance(q, Mul)):
        ops.append(q)
    else:
        for arg in q.args:
            if is_ladder_contained(arg):
                ops.append(arg)
            else:
                scalars.append(arg)
                
    return Mul(*scalars), Mul(*ops) 
            # Empty sequence is evaluated to 1.
                
def expval_preprocess(val):
    val = sympify(val).expand().factor()

    common_factor = 1
    if isinstance(val, Mul) and isinstance(val.args[1], Add):
        common_factor = val.args[0]
        assert not(isinstance(common_factor, Add))
        
        val = Mul(*val.args[1:]).expand()
    
    scalars = []
    ops = []
    
    if not(isinstance(val, Add)):
        scalar, op = get_scalars_and_ops(val)
        scalars.append(scalar)
        ops.append([op])
    else:
        for arg in val.args:
            scalar, op = get_scalars_and_ops(arg)
            scalars.append(scalar)
            ops.append([op])   
                        # This list is appended to show
                        # multiplication of expectation
                        # values, defined in the main 
                        # class.
    
    return common_factor, scalars, ops

def sympify_and_expval(q):
    if not(isinstance(q, expval)):
        return expval(sympify(q))
    else:
        return q

class expval(Expr):
    """
    Expectation value object.
    """
    
    def __init__(self, val, normal_order = False):
        super().__init__()
        val = sympify(val)
        if normal_order:
            val = normal_order(val)
        self.common_factor, \
            self.scalars, \
            self.ops = expval_preprocess(val)

    def _latex(self, printer):
                        
        out = "" 
        for i in range(len(self.scalars)):                            
            
            if self.scalars[i] not in [1, Integer(1)]:
                out += r"%s"%(printer.doprint(self.scalars[i]))
            
            j = 0
            while True:
                # self.ops is a list of summands, where each element
                # is a list of factors of expectation values. Mul
                # cannot be used here as we want to raise the power
                # of the expval object instead of what it contains.
                print(i, j, self.ops)
                to_write = self.ops[i][j]
                
                if to_write in [1, Integer(1)]:
                    if len(self.ops[i]) == 1:
                        out += printer.doprint(Integer(1))
                        break
                    j += 1
                    continue
                                    
                pow = 1
                if j != (len(self.ops[i])-1):
                # The last quantity that is not one is 
                    for q in self.ops[i][j+1:]:
                        if q != to_write:
                            break
                        pow += 1
                
                out += r"\left\langle %s \right\rangle %s"%(printer.doprint(to_write),
                                                            
                                                            r"^{%s}"%(printer.doprint(pow)) 
                                                                if pow>1 
                                                                else r"")
                
                j += pow
                # Jump to the factor next to the last q that is the same as to_write.
                
                if j == len(self.ops[i]):
                    break
                    
            out += r"+"     
    
        # Remove unnecessary plus signs.
        i = 0
        while True:
            if i == (len(out)-1):
                break
            if out[i] != "+":
                i += 1
                continue
            if out[i+1] == "-":
                out = out[:i] + out[i+1:]
        out = out[:-1]
            
        if self.common_factor != 1:
            out = r"%s\left["%(printer.doprint(self.common_factor)) \
                    + out + r"\right]"    
    
        return out
        
    def __add__(self, x):
        x = sympify_and_expval(x)
        
        share_common_factor = self.common_factor == x.common_factor
        res_common_factor = self.common_factor if share_common_factor \
                                                else 1
        if not(share_common_factor):
            self.scalars = [q*self.common_factor for q in self.scalars]
            x_scalars = [q*x.common_factor for q in x_scalars]
        
        res_scalars = []
        res_ops = []
        
        i_not_added = list(range(len(self.scalars)))
        j_not_added = list(range(len(x.scalars)))
        
        for i, val_i in enumerate(self.scalars):        
            for j, val_j in enumerate(x.scalars):
                if self.ops[i] == x.ops[j]:
                    res_scalars.append(val_i+val_j)
                    res_ops.append(self.ops[i])
                    i_not_added.pop(i)
                    j_not_added.pop(j)
        
        for i in i_not_added:
            res_scalars.append(self.scalars[i])
            res_ops.append(self.ops[i])
        for j in j_not_added:
            res_scalars.append(self.scalars[j])
            res_ops.append(self.ops[j])
            
        res =  expval(1)
        res.common_factor = res_common_factor
        res.scalars = res_scalars
        res.ops = res_ops
        
        return res
            
    def __mul__(self, x):
        x = sympify_and_expval(x)
        
        self_scalars = [self.common_factor*q for q in self.scalars]
        x_scalars = [x.common_factor*q for q in x.scalars]
        
        res_scalars = [a*b for a in self_scalars for b in x_scalars]
        res_common_factor = gcd(*res_scalars) if len(res_scalars)>1 else 1
        res_scalars = [q/res_common_factor for q in res_scalars]
        
        res_ops = []
        for op_a in self.ops:
            for op_b in x.ops:
                res_ops.append(op_a+op_b)
        
        res = expval(1)
        res.common_factor = res_common_factor
        res.scalars = res_scalars
        res.ops = res_ops
        
        return res
    
    def __pow__(self, x):
        res = 1
        for _ in range(x):
            res *= self
            
        return res 