from __future__ import (absolute_import, division, print_function, unicode_literals)

import starlink.Ast as Ast

from .region import ASTRegion
from ..mapping import ASTMapping
from ..mapping import ASTFrame

__all__ = ["ASTBox"]

CENTER_CORNER = 0
CORNER_CORNER = 1

class ASTBox(ASTRegion):
	'''
	ASTBox is an ASTRegion that represents a box with sides parallel to the axes of an ASTFrame.
	
	There are two accepted signatures for creating an ASTBox:
	
	b = ASTBox(frame, cornerPoint, cornerPoint2)
	b = ASTBox(frame, cornerPoint, centerPoint)
	
	Points can be any two element container, e.g. (1,2), [1,2], np.array([1,2])
	
	The 'frame' parameter may either be an ASTFrame object or a starlink.Ast.frame object.
	
	A Box is similar to an Interval, the only real difference being that the Interval
	class allows some axis limits to be unspecified. Note, a Box will only look like a box
	if the Frame geometry is approximately flat. For instance, a Box centered close to a pole
	in a SkyFrame will look more like a fan than a box (the Polygon class can be used to
	create a box-like region close to a pole).
	
	self.astObject is of type starlink.Ast.Box.
	'''
	def __init__(self, frame=None, cornerPoint=None, cornerPoint2=None, centerPoint=None):
		#self.astFrame = frame
		self._uncertainty = 4.848e-6 # defaults to 1 arcsec
		#self._ast_box = None
		
		# input forms:
		#    0: box specified by center point and any corner point
		#    1: box specified by a corner and its oppsite corner
		input_form = None
		
		# check valid combination of parameters
		# -------------------------------------
		if frame is None:
			raise Exception("A frame must be specified when creating an ASTBox.")
		else:
			if isinstance(frame, ASTFrame):
				self.frame = frame
			elif isinstance(frame, starlink.Ast.frame):
				self.frame = ASTFrame(frame=frame)

		if all([cornerPoint,centerPoint]) or all([cornerPoint,cornerPoint2]):
			if centerPoint is None:
				input_form = CORNER_CORNER
			else:
				input_form = CENTER_CORNER
		else:
			raise Exception("Either 'cornerPoint' and 'centerPoint' or 'cornerPoint' " + \
							"and 'cornerPoint2' must be specified when creating an ASTBox.")
		
		if input_form == CENTER_CORNER:
			p1 = [centerPoint[0], centerPoint[1]]
			p2 = [cornerPoint[0], cornerPoint[1]]
		else:
			p1 = [cornerPoint[0], cornerPoint[1]]
			p2 = [cornerPoint2[0], cornerPoint2[1]]
		
		# Box args: :frame,form,point1,point2,unc=None,options=None  <-- note which are keyword args & which not
		# AstBox( starlink.Ast.Frame(2), [0,1], 
		self.astObject = Ast.Box( self.frame.astObject, input_form, p1, p2, unc=self.uncertainty )

	@property
	def uncertainty(self):
		'''
		Uncertainties associated with the boundary of the Box.
					
		The uncertainty in any point on the boundary of the Box is found by
		shifting the supplied "uncertainty" Region so that it is centered at
		the boundary point being considered.  The area covered by the shifted
		uncertainty Region then represents the uncertainty in the boundary
		position.  The uncertainty is assumed to be the same for all points.
		'''
		return self._uncertainty
			
	@uncertainty.setter
	def uncertainty(self, unc):
		raise Exception("Setting the uncertainty currently doesn't do anything.")
		self._uncertainty = unc
	
	def mapRegionMesh(self, mapping=None, frame=None):
		'''
		Returns a new ASTRegion that is the same as this one but with the specified coordinate system.
		
		Parameters
		----------
		mapping : `~cornish.mapping.ASTMapping` class
			The mapping to transform positions from the current ASTRegion to those specified by the given frame.
		frame : `~cornish.frame.ASTFrame` class
			Coordinate system frame to convert the current ASTRegion to.
			
		Returns
		-------
		region : `ASTRegion`
			A new region object covering the same area but in the frame specified in `frame`.
			
		Raises
		------
		Exception
			An exception is raised for missing parameters.
		'''
		if mapping is None or frame is None:
			raise Exception("A mapping and frame is required to be passed to 'mapRegion'.")

		# check it's the correct type
		if not isinstance(mapping, ASTMapping):
			raise Exception("The object passed to 'mapping' needs to be an ASTMapping object.")
		
		if not isinstance(frame, ASTFrame):
			raise Exception("The object passed to 'frame' needs to be an ASTFrame object.")
		
		self.astObject.mapregionmesh( mapping, frame )

			
			
		