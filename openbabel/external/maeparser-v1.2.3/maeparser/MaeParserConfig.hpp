#pragma once

#ifndef STATIC_MAEPARSER

#ifdef IN_MAEPARSER

#ifdef WIN32
#define EXPORT_MAEPARSER __declspec(dllexport)
#else
#define EXPORT_MAEPARSER __attribute__((visibility("default")))
#endif

#else

#ifdef WIN32
#define EXPORT_MAEPARSER __declspec(dllimport)
#else
#define EXPORT_MAEPARSER
#endif

#endif

#else

#define EXPORT_MAEPARSER

#endif
