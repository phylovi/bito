// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"


GPEngine::GPEngine(SitePattern site_pattern, size_t pcss_count)
    : site_pattern_(std::move(site_pattern)), pcss_count_(pcss_count) {}

