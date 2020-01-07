#include "lib.hpp"

int main() {
    HDRHist hist;
    hist.add_value(1000000);
    auto ccdf = hist.ccdf();
    for (
        std::optional<CcdfElement> ccdf_el = ccdf.next();
        ccdf_el.has_value();
        ccdf_el = ccdf.next()) {

    }
}
