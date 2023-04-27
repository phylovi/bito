#pragma once

#include "../src/driver.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Driver") {
  Driver driver;

  std::vector<std::string> newicks = {
      "(a:0,b:0,c:0,d:0):0;",
      "((b:0,a:0):0,c:0):0;",
      "((a:1.1,b:2):0.4,c:3):0;",
      "(x:0,(a:1.1,(b:2,(quack:0.1,duck:0):0):0):0,c:3):1.1;",
  };
  for (const auto& newick : newicks) {
    auto collection = driver.ParseString(newick);
    CHECK_EQ(newick, collection.Trees()[0].Newick(collection.TagTaxonMap()));
  }
  driver.Clear();
  // Note that the order of the taxa is given by the order in the translate table, not
  // by the short names. We use that here to make sure that the ordering of the taxa is
  // the same as that in the newick file below so that they can be compared.
  auto nexus_collection = driver.ParseNexusFile("data/DS1.subsampled_10.t.reordered");
  CHECK_EQ(nexus_collection.TreeCount(), 10);
  driver.Clear();
  auto newick_collection = driver.ParseNewickFile("data/DS1.subsampled_10.t.nwk");
  CHECK_EQ(nexus_collection, newick_collection);
  driver.Clear();
  auto newick_collection_gz =
      driver.ParseNewickFileGZ("data/DS1.subsampled_10.t.nwk.gz");
  CHECK_EQ(nexus_collection, newick_collection_gz);
  driver.Clear();
  auto five_taxon = driver.ParseNewickFile("data/five_taxon_unrooted.nwk");
  std::vector<std::string> correct_five_taxon_names({"x0", "x1", "x2", "x3", "x4"});
  CHECK_EQ(five_taxon.TaxonNames(), correct_five_taxon_names);
  // Check that we can parse BEAST trees with [&comments], and that the different
  // formatting of the translate block doesn't trip us up.
  auto beast_nexus = driver.ParseNexusFile("data/test_beast_tree_parsing.nexus");
  // These are the taxa, in order, taken directly from the nexus file:
  StringVector beast_taxa = {
      "aDuckA_1976",     "aDuckB_1977",   "aItaly_1987",   "aMallard_1985",
      "hCHR_1983",       "hCambr_1939",   "hFortMon_1947", "hKiev_1979",
      "hLenin_1954",     "hMongol_1985",  "hMongol_1991",  "hNWS_1933",
      "hPR_1934",        "hSCar_1918.00", "hScot_1994",    "hSuita_1989",
      "hUSSR_1977",      "sEhime_1980",   "sIllino_1963",  "sIowa_1930",
      "sNebrask_1992",   "sNewJers_1976", "sStHya_1991",   "sWiscons_1961",
      "sWiscons_1.998e3"};
  CHECK_EQ(beast_nexus.TaxonNames(), beast_taxa);
  // Check that we got the whole tree.
  for (const auto& [topology, count] : beast_nexus.TopologyCounter()) {
    std::ignore = count;
    CHECK_EQ(topology->LeafCount(), beast_taxa.size());
  }
  auto beast_nexus_gz =
      driver.ParseNexusFileGZ("data/test_beast_tree_parsing.nexus.gz");
  CHECK_EQ(beast_nexus, beast_nexus_gz);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
