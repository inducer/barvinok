AC_DEFUN([AX_UNORDERED_MAP], [
	AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	AC_TRY_COMPILE([#include <unordered_map>], [using std::unordered_map;],
		AC_DEFINE([HAVE_UNORDERED_MAP], [],
			[define if std::unordered_map is available]))
	AC_LANG_RESTORE
])
