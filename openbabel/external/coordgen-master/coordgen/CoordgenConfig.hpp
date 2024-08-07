#pragma once

#ifndef STATIC_COORDGEN

#ifdef IN_COORDGEN
#ifdef WIN32
#define EXPORT_COORDGEN __declspec(dllexport)
#else
#define EXPORT_COORDGEN __attribute__((visibility("default")))
#endif

#else

#ifdef WIN32
#define EXPORT_COORDGEN __declspec(dllimport)
#else
#define EXPORT_COORDGEN
#endif

#endif

#else

#define EXPORT_COORDGEN

#endif
