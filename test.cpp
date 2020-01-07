#include "hdrhist.hpp"
#include <iostream>

int main() {
    HDRHist hist;
    hist.add_value(1000000);
    hist.add_value(2000000);
    hist.add_value(3000000);
    auto ccdf = hist.ccdf();
    for (
        std::optional<CcdfElement> ccdf_el = ccdf.next();
        ccdf_el.has_value();
        ccdf_el = ccdf.next()) {

        std::cerr << ccdf_el->value << " " << ccdf_el->fraction << std::endl;
    }

    auto quantiles_f = std::vector { .25f, .5f, .75f };
    auto quantiles = hist.quantiles(quantiles_f);
    for (
        auto quantile_el = quantiles.next();
        quantile_el.has_value();
        quantile_el = quantiles.next()) {

        std::cerr << quantile_el->quantile << " " << quantile_el->lower_bound << " " << quantile_el->upper_bound << std::endl;
    }
}
