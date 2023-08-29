#!/usr/bin/env python3      

import re
from sympy import Matrix, lcm, Rational
elementList=[]
elementMatrix=[]
print("please input your reactants, this is case sensitive")
print("your input should look like: H2O+Ag3(Fe3O)4")
# reactants=input("Reactants: ")
print("please input your products, this is case sensitive")
# products=input("Products: ")
# reactants=reactants.replace(' ', '').split("+")
# products=products.replace(' ', '').split("+")
reactants= ["LaH3", "CeH2", "ThH3", "YH3", "BeH2", "H"]
products=["(LaCe)3ThYBe4H32"]
def addToMatrix(element, index, count, side):
    print(f"index={index}, elementMatrix={elementMatrix}"); input()
    if(index == len(elementMatrix)):
       elementMatrix.append([])
       for x in elementList:
            elementMatrix[index].append(0)
    print(f"第一次：elementMatrix={elementMatrix}")
    print(f"element={element}")
    if(element not in elementList):
        elementList.append(element)
        for i in range(len(elementMatrix)):
            elementMatrix[i].append(0)
    column=elementList.index(element)
    elementMatrix[index][column]+=count*side
    print(f"第二次：elementMatrix={elementMatrix}")

def findElements(segment, index, multiplier, side):
    print(f"segment={segment}")
    elementsAndNumbers=re.split('([A-Z][a-z]?)',segment)
    print(elementsAndNumbers, "elementsAndNumbers"); input()
    i=0
    while(i<len(elementsAndNumbers)-1):#last element always blank
          i+=1
          if(len(elementsAndNumbers[i])>0):
            print(f"i={i}, elementsAndNumbers[i+1].isdigit()={elementsAndNumbers[i+1].isdigit()}")
            if(elementsAndNumbers[i+1].isdigit()):
                count=int(elementsAndNumbers[i+1])*multiplier
                addToMatrix(elementsAndNumbers[i], index, count, side)
                i+=1
            else:
                addToMatrix(elementsAndNumbers[i], index, multiplier, side)        
    
def compoundDecipher(compound, index, side):
    segments=re.split('(\([A-Za-z0-9]*\)[0-9]*)',compound)
    input(segments)
    for segment in segments:
        if segment.startswith("("):
            segment=re.split('\)([0-9]*)',segment) # 这部分用圆括号括起来，表示一个捕获组。它匹配零个或多个数字（0-9）
            print(segment)
            multiplier=int(segment[1])
            segment=segment[0][1:]
        else:
            multiplier=1
        findElements(segment, index, multiplier, side)


if __name__ == "__main__":
    for i in range(len(reactants)):
        compoundDecipher(reactants[i],i,1)
    for i in range(len(products)):
        compoundDecipher(products[i],i+len(reactants),-1)
    elementMatrix = Matrix(elementMatrix)
    print(f"elementMatrix={elementMatrix}")
    elementMatrix = elementMatrix.transpose()
    print(f"aftertranpose-elementMatrix={elementMatrix}")
    solution=elementMatrix.nullspace()[0]
    print(f"elementMatrix.nullspace()={elementMatrix.nullspace()}")
    print(f"solution={solution}")
    multiple = lcm([val.q for val in solution]) # 获得所有系数的分母的值，然后取它们的最小公倍数
    # 在 SymPy（一个用于符号数学计算的Python库）中，.q 是用来获取有理数对象（Rational）的分母的属性。
    # 有理数是一个分数形式的数值，例如 3/4，其中分子为3，分母为4。
    # 通过 .q 属性，你可以获取有理数对象的分母部分。
    solution = multiple*solution

    coEffi=solution.tolist()
    print(f"coEffi={coEffi}"); input()
    output=""
    for i in range(len(reactants)):
        output+=str(coEffi[i][0])+reactants[i]
        if i<len(reactants)-1:
            output+=" + "
    output+=" -> "
    for i in range(len(products)):
        output+=str(coEffi[i+len(reactants)][0])+products[i]
        if i<len(products)-1:
            output+=" + "
    print(output)

    

