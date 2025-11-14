#pragma once

#include <vector>
#include <string>

// floating point type
using FloatType = double;

// index type
using IndexType = size_t;

// type to be used for a scalar field
using FieldType = typename std::vector<std::vector<FloatType>>;

// enum and name (identifiers) for primitive variables
enum PV { U = 0, V, P, U_OLD, V_OLD, P_OLD, U_PICARD_OLD, V_PICARD_OLD, P_PICARD_OLD };
const static std::vector<std::string> pvNames = { "u", "v", "p", "uOld", "vOld", "pOld", "uPicardOld", "vPicardOld", "pPicardOld" };