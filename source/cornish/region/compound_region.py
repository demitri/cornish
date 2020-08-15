
from typing import Iterable, Union

import starlink.Ast as Ast

from .region import ASTRegion

__all__ = ["ASTCompoundRegion"]

class ASTCompoundRegion(ASTRegion):
	'''
	A region that is created as the composite of multiple regions.
	
	Regions are composited two at a time in the order they are supplied,
	e.g. if regions=[r1, r2, r3, r4]
	the result would be
	
	region = compound( compound( compound(r1, r2), r3 ) r4 )
	
	all using the same operator as specified.
	
	:param regions: a list of regions to compound
	:param operation: one of `starlink.Ast.AND, `starlink.Ast.OR`, `starlink.Ast.XOR`
	'''
	def __init__(self, ast_object=None, regions:Iterable[Union[ASTRegion, Ast.Region]]=None, operation:int=None):
		if ast_object:
			if any([regions, operation]):
				raise ValueError("If the ast_object is specified, no other parameter is accepted.")
			else:
				super().__init__(ast_object=ast_object)
				return
				
		if regions is None or len(regions) < 2:
			raise ValueError("The 'regions' parameter must contain at least two regions.")
		for r in regions:
			if not isinstance(r, (ASTRegion, Ast.Region)):
				raise ValueError("The regions provided must be of type ASTRegion or starlink.Ast.Region.")
		
		r1 = None
		r2 = None
		compound_region = None
		
		while len(regions) > 0:
			if r1 is None:
				r1 = regions.pop(0) # get first item
				r2 = regions.pop(0)
			else:
				r1 = compound_region
				r2 = regions.pop(0)
				
			if isinstance(r1, ASTRegion):
				r1 = r1.astObject
			if isinstance(r2, ASTRegion):
				r2 = r2.astObject

			compound_region = Ast.CmpRegion( r1, r2, oper=operation ) # todo: 'series' parameter?
			if compound_region is None:
				# an error occurred
				return None
		
		super().__init__(ast_object=compound_region)
