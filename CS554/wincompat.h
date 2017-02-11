#include <stdio.h>
#include <tchar.h>

/** Code borrowed from http://www.wincli.com/?p=72 */

#ifndef __WINCOMPAT_H__

	#define __WINCOMPAT_H__

	int wstrlen(_TCHAR * wstr);
	char * wstrdup(_TCHAR * wSrc);
	char ** allocate_argn (int argc, _TCHAR* argv[]);
	void release_argn(int argc, char ** nargv);

#endif
