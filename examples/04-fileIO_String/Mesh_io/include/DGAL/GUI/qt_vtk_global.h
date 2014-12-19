#ifndef QT_VTK_GLOBAL_H
#define QT_VTK_GLOBAL_H

#include <Qt/qglobal.h>
#include <QtGui>

#ifdef QT_VTK_LIB
# define QT_VTK_EXPORT Q_DECL_EXPORT
#else
# define QT_VTK_EXPORT Q_DECL_IMPORT
#endif

#endif // QT_VTK_GLOBAL_H
