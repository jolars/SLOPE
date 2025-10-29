/**
 * @file
 * @brief Main header for the SLOPE library
 *
 * This header provides the main public interface for the SLOPE (Sorted L-One
 * Penalized Estimation) library.
 */
#pragma once

// Core SLOPE functionality
#include <slope/slope.h>
#include <slope/slope_fit.h>
#include <slope/slope_path.h>

// Cross-validation
#include <slope/cv.h>

// Proximal operator
#include <slope/sorted_l1_norm.h>
