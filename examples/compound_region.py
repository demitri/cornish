
import numpy as np
import starlink.Ast as Ast
from cornish import ASTCircle, ASTICRSFrame, ASTCompoundRegion

frame = ASTICRSFrame()
c1 = ASTCircle(frame=frame, center=[12, 0], radius=12.)
c2 = ASTCircle(frame=frame, center=[5, 5], radius=12.)

cr = ASTCompoundRegion(regions=[c1,c2], operation=Ast.AND)
print(cr)
