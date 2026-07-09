
'''
Small source/sink adapter objects for AST channel I/O.

AST channels read and write text through objects that provide an
``astsource()`` method (return the next input line, or None when done) and an
``astsink(line)`` method (receive one output line). These two adapters cover
the common case of reading from / writing to Python lists of strings so that
callers never need to know the protocol exists.
'''

__all__ = ['ListSink', 'ListSource']

class ListSink:
	''' Collects lines written by an AST channel into ``self.lines``. '''
	def __init__(self):
		self.lines = list()

	def astsink(self, line:str):
		self.lines.append(line)

class ListSource:
	''' Feeds lines of text (e.g. a serialized AST object) to an AST channel. '''
	def __init__(self, lines):
		if isinstance(lines, str):
			lines = lines.splitlines()
		self.lines = list(lines)
		self._index = 0

	def astsource(self):
		if self._index < len(self.lines):
			line = self.lines[self._index]
			self._index += 1
			return line
		return None
