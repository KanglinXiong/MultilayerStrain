#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################


import copy
import Misc


class ArrayOp:

    '''
        basic operation for matrix and vector
        ranks are not checked
    '''

    def __init__(self, tolerance = 1e-10):
        self.tolerance = tolerance
        pass

    @staticmethod
    def vecDotVec(v1, v2):
        indices = range(len(v1))
        return sum(map(lambda i: v1[i]*v2[i], indices)) 
        pass

    @staticmethod
    def vecAddVec(v1, v2):
        indices = range(len(v1))
        return list(map(lambda i: v1[i]+v2[i], indices)) 
        pass

    @staticmethod
    def vecMltSca(v, s):
        indices = range(len(v))
        return list(map(lambda i: v[i]*s, indices)) 
        pass
 
    @classmethod
    def matMltSca(cls, m, s):
        indices = range(len(m))
        return list(map(lambda i: cls.vecMltSca(m[i], s), indices))
        pass
 
    @classmethod
    def matAddMat(cls, m1, m2):
        indices = range(len(m1))
        return list(map(lambda i: cls.vecAddVec(m1[i], m2[i]), indices))
        pass

    @classmethod
    def matDotVec(cls, m, v):
        indices = range(len(m))
        return list(map(lambda i: cls.vecDotVec(m[i], v), indices))
        pass

    @staticmethod
    def matCol(m, j):
        indices = range(len(m))
        return list(map(lambda i: m[i][j], indices))
        pass

    @classmethod
    def matTranspose(cls, m):
        indices = range(len(m[0]))
        return list(map(lambda j: cls.matCol(m, j), indices))
        pass

    @classmethod
    def matDotMat(cls, m1, m2):
        if(len(m1) != len(m2[0])):
            raise Exception("Invalid matrix rank!")
        m2T = cls.matTranspose(m2)
        indices = range(len(m2T))
        return cls.matTranspose(list(map(lambda i: cls.matDotVec(m1, m2T[i]), indices)))
        pass

    @classmethod
    def identityMat(cls, n):
        indices = range(n)
        m = list(map(lambda i: range(n), indices))
        m = cls.matMltSca(m, 0.0)
        for i in indices:
            m[i][i] = 1.0
        return m
        pass

    @classmethod
    def diagMat(cls, v):
        m = cls.identityMat(len(v))
        indices = range(len(v))
        for i in indices:
            m[i][i] = v[i]
        return m
        pass


class LinearEq:

    '''
        Solve linear equation m.x == b
        can also find inverse of m
    '''
   
    def __init__(self, m, b = [], matInverseFlag = True, relativeTolerance = 1e-8):
        if(len(m) != len(m[0])):
            raise Exception("Invalid matrix!")
        self.matrix = copy.deepcopy(m)
        self.rank = len(m)
        if(len(b) != len(m)):
            self.vector =  ArrayOp.vecMltSca(list(range(self.rank), 0.0))
        else:
            self.vector = copy.deepcopy(b)
        self.matInverseFlag = matInverseFlag
        maxAbsMatEle = max(map(lambda i: (abs(max(m[i]))+abs(min(m[i]))), range(self.rank)))
        self.absTolerance = abs(relativeTolerance * maxAbsMatEle)
        self.workMat = list(map(lambda i: [], range(self.rank)))
        self.concatenate()
        self.invMat = None
        self.root = None
        pass

    def concatenate(self):
        ''' laterally concatenate m, b, identity ''' 
        if self.matInverseFlag:
            idMat = ArrayOp.identityMat(self.rank)
        for i in range(self.rank):
            self.workMat[i].extend(self.matrix[i])
            self.workMat[i].extend([self.vector[i]])
            if(self.matInverseFlag):
                self.workMat[i].extend(idMat[i])
        pass

    @staticmethod
    def maxAbsElePos(vec):
        idx = 0
        for i in range(len(vec)):
            if(abs(vec[idx]) < abs(vec[i])):
                idx = i
        return idx

    def precondition(self):
        ''' before solving, deal with zero diagonal elements '''
        for i in range(self.rank):
            if(abs(self.workMat[i][i]) > self.absTolerance): continue
            # find a col to overwrite the zero
            colVec = ArrayOp.matCol(self.workMat, i) 
            rowNum = self.maxAbsElePos(colVec)
            self.workMat[i] = ArrayOp.vecAddVec(self.workMat[i], self.workMat[rowNum])
        pass

    def runCol(self, colNum):
        ''' make diagonal element 1 and other 0, by scaling and row operations '''
        if(colNum < 0 or colNum >= self.rank):
            raise Exception("Wrong column number!")
        # scale diagonal element to 1
        scaleFactor = 1.0/self.workMat[colNum][colNum]
        self.workMat[colNum] = ArrayOp.vecMltSca(self.workMat[colNum], scaleFactor)
        # make other element 0
        for i in range(self.rank):
            if i==colNum: continue
            scaleFactor = 0.0 - self.workMat[i][colNum]
            self.workMat[i] = ArrayOp.vecAddVec(\
                self.workMat[i],\
                ArrayOp.vecMltSca(self.workMat[colNum], scaleFactor))
        pass

    def solve(self):
        ''' column by column '''
        self.precondition()
        for i in range(self.rank):
            self.runCol(i)
        self.invMat = list(map(\
            lambda i: self.workMat[i][0+self.rank+1 : 0+self.rank+1+self.rank],\
            range(self.rank)))
        self.root = ArrayOp.matCol(self.workMat,0+self.rank)
        pass

    def error(self):
        lhs = ArrayOp.matDotVec(self.matrix, self.root)
        rhs = ArrayOp.vecMltSca(self.vector, -1.0)
        vec = ArrayOp.vecAddVec(lhs, rhs)
        return sum(map(lambda x: abs(x), vec))/len(vec)

    def getRoot(self):
        return self.root

    def getInverseMatrix(self):
        return self.InvMat

if __name__ == "__main__":
    print("test ArrayOp")
    op = ArrayOp 
    v1 = [1.0, 2]
    v2 = [3, 4.0]
    print("v1 and v2:")
    print(v1)
    print(v2)
    print("v1 dot v2")
    print(op.vecDotVec(v1, v2))
    print("v1 add v2")
    print(op.vecAddVec(v1, v2))
    print("[v1, v2]*10")
    Misc.display(op.matMltSca([v1, v2], 10))
    print("[v1, v1] + [v2, v2]")
    Misc.display(op.matAddMat([v1, v1], [v2, v2]))
    print("[v1, v2].v2")
    print(op.matDotVec([v1, v2], v2))
    print("transpose [v1, v2, v1]")
    Misc.display(op.matTranspose([v1, v2, v1]))
    print("[v1, v2].[v2, v1]")
    Misc.display(op.matDotMat([v1, v2], [v2, v1]))
    print("identityMat(4)")
    Misc.display(op.identityMat(4))
    print("diagMat(range(1,11))")
    Misc.display(op.diagMat(list(range(1,11))))
    print("range(1,11).range(2,12)")
    Misc.display(op.matDotMat(op.diagMat(list(range(1,11))), op.diagMat(list(range(2,12)))))

if __name__ == "__main__":
    print("\ntest LinearEq")
    m = [[0, 10], [3, 8]]
    b = [9, 5]
    eq = LinearEq(m, b, True)
    print("abs tolerance")
    print(eq.absTolerance)
    print("inital m, b")
    Misc.display(m)
    print(b)
    eq.solve()
    print("post run m, b, I")
    Misc.display(eq.workMat)
    print("invMat, root")
    Misc.display(eq.invMat)
    print(eq.root)
    Misc.display("m.invMat")
    Misc.display(ArrayOp.matDotMat(m, eq.invMat))
    print("m.root")
    print(ArrayOp.matDotVec(m, eq.root))
    print(eq.error())


