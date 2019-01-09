/**********************************************************************
invector.h - Convenience function to determine std::vector membership
***********************************************************************/

#ifndef INVECTOR_H
#define INVECTOR_H

#include <vector>
#include <algorithm>

// Keeping this inline template in the .h file to avoid linker complications.
// For more details, see https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
template<typename T>
bool inVector(const T &element, const std::vector<T> &vec) {
	// Test if an element is a given member of a vector
	if (std::find(vec.begin(), vec.end(), element) != vec.end()) {  // idiomatic C++
		return true;
	} else {
		return false;
	}
}

#endif // INVECTOR_H

//! \file invector.h
//! \brief Convenience function to determine std::vector membership
