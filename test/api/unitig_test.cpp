// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "../app_test.hpp"
#include "minimiser.hpp"
#include "misc/fill_hash_table.hpp"
#include "misc/stream.hpp"

struct unitig_test : public app_test
{
    minimiser_arguments args{};
    robin_hood::unordered_node_map<uint64_t, uint16_t> hashes{};
    robin_hood::unordered_set<uint64_t> include{}, exclude{};

    unitig_test()
    {
        args.k = 4;
        args.shape = seqan3::ungapped{4};
        args.w_size = seqan3::window_size{4};
        args.s = seqan3::seed{0};
    }

    void read(std::string const & fasta, uint8_t cutoff = 0, bool only_include = false)
    {
        std::istringstream input{fasta};
        sequence_file_with_id_t fin{input, seqan3::format_fasta{}};
        fill_hash_table_unitigs(args, fin, hashes, include, exclude, only_include, cutoff);
    }
};

TEST_F(unitig_test, header_abundance_replaces_occurrence_counts)
{
    read(">SRR1_0 LN:i:8 ka:f:3.4 L:+:1:+\nAAAAAAAA\n"
         ">SRR1_1 km:f:7.5\nTTTTTTTT\n"
         ">SRR1_2 ka:f:2\nAAAA\n"
         ">SRR1_3 ka:f:70000\nCCCC\n"
         ">SRR1_4 ka:f:0\nACGT\n"
         ">SRR1_5 ka:f:12\nACG\n");
    ASSERT_EQ(hashes.size(), 2u);
    EXPECT_EQ(hashes.at(0), 8);      // AAAA and its reverse complement share a hash.
    EXPECT_EQ(hashes.at(85), 65534); // CCCC; saturate before conversion to uint16_t.
}

TEST_F(unitig_test, cutoff_and_filters)
{
    std::string const fasta{">a ka:f:3\nAAAA\n>b ka:f:4\nCCCC\n>c ka:f:5\nACGT\n"};
    read(fasta, 3);
    EXPECT_EQ(hashes.size(), 2u);
    EXPECT_FALSE(hashes.contains(0));
    EXPECT_EQ(hashes.at(85), 4);
    EXPECT_EQ(hashes.at(27), 5);

    hashes.clear();
    exclude.insert(85);
    read(fasta, 3);
    ASSERT_EQ(hashes.size(), 1u);
    EXPECT_EQ(hashes.at(27), 5);

    hashes.clear();
    read(fasta, 0, true); // An empty include set includes nothing.
    EXPECT_TRUE(hashes.empty());
    include.insert(85);
    read(fasta, 0, true);
    ASSERT_EQ(hashes.size(), 1u);
    EXPECT_EQ(hashes.at(85), 4);
}

TEST_F(unitig_test, fractional_abundances)
{
    read(">a ka:f:0.49\nAAAA\n>b ka:f:0.5\nCCCC\n>c ka:f:3.4\nACGT\n");
    ASSERT_EQ(hashes.size(), 2u);
    EXPECT_FALSE(hashes.contains(0));
    EXPECT_EQ(hashes.at(85), 1);
    EXPECT_EQ(hashes.at(27), 3);
}

TEST_F(unitig_test, invalid_abundance)
{
    for (std::string const field :
         {"", "KC:i:4", "ka:f:", "ka:f:abc", "ka:f:2junk", "ka:f:-1", "ka:f:nan", "ka:f:inf", "ka:f:1e999", "km:f:-2"})
    {
        SCOPED_TRACE(field);
        EXPECT_THROW(read(">SRR1_0 " + field + "\nAAAA\n"), std::invalid_argument);
    }
}

TEST_F(unitig_test, grouped_files_and_ram_mode)
{
    std::ofstream{"a.fa"} << ">SRR1_0 ka:f:3\nAAAA\n>SRR1_1 ka:f:9\nCCCC\n";
    std::ofstream{"b.fa"} << ">SRR1_2 km:f:7\nTTTT\n";
    args.threads = 2;
    minimiser_file_input_arguments input_args{};
    input_args.unitigs = true;
    input_args.samples = {2};
    for (bool const ram : {false, true})
    {
        input_args.ram_friendly = ram;
        std::vector<uint8_t> cutoffs{};
        minimiser({"a.fa", "b.fa"}, args, input_args, cutoffs);
        robin_hood::unordered_node_map<uint64_t, uint16_t> result;
        read_binary("a.minimiser", result);
        ASSERT_EQ(result.size(), 2u);
        EXPECT_EQ(result.at(0), 7);
        EXPECT_EQ(result.at(85), 9);
        ASSERT_EQ(cutoffs.size(), 1u);
        EXPECT_EQ(cutoffs[0], 0);

        estimate_ibf_arguments stored_args{};
        uint64_t size{};
        uint8_t cutoff{255};
        read_binary_start(stored_args, "a.minimiser", size, cutoff);
        EXPECT_EQ(size, 2u);
        EXPECT_EQ(cutoff, 0);
        EXPECT_EQ(stored_args.k, 4);
        EXPECT_EQ(stored_args.w_size.get(), 4u);
    }

    // Also exercise parallel processing of independent samples.
    input_args.samples.clear();
    input_args.ram_friendly = false;
    std::vector<uint8_t> cutoffs{};
    minimiser({"a.fa", "b.fa"}, args, input_args, cutoffs);
    read_binary("a.minimiser", hashes);
    EXPECT_EQ(hashes.at(0), 3);
    hashes.clear();
    read_binary("b.minimiser", hashes);
    ASSERT_EQ(hashes.size(), 1u);
    EXPECT_EQ(hashes.at(0), 7);
}

TEST_F(unitig_test, explicit_cutoff)
{
    std::ofstream{"a.fa"} << ">a ka:f:3\nAAAA\n>b ka:f:4\nCCCC\n";
    minimiser_file_input_arguments input_args{};
    input_args.unitigs = true;
    std::vector<uint8_t> cutoffs{3};
    minimiser({"a.fa"}, args, input_args, cutoffs);
    read_binary("a.minimiser", hashes);
    ASSERT_EQ(hashes.size(), 1u);
    EXPECT_EQ(hashes.at(85), 4);
}

TEST_F(unitig_test, cli_build_index_from_logan_unitigs)
{
    std::ofstream{"logan.fa"} << ">SRR1_0 ka:f:7.5\n"
                              << std::string(40, 'A') << "C\n>SRR1_1 km:f:12\n"
                              << std::string(31, 'C') << "A\n";
    app_test_result result = execute_app("minimiser --unitigs -k 31 -w 31 --seed 0 logan.fa");
    ASSERT_EQ(result.exit_code, 0) << result.err;
    read_binary("logan.minimiser", hashes);
    ASSERT_EQ(hashes.size(), 4u);
    EXPECT_EQ(hashes.at(0), 8);

    result = execute_app("ibfmin logan.minimiser -e 5 -e 10 -f 0.0001 -o index_");
    ASSERT_EQ(result.exit_code, 0) << result.err;
    seqan3::interleaved_bloom_filter<> index;
    load_ibf(index, "index_IBF_5");
    auto low = index.membership_agent();
    EXPECT_TRUE(low.bulk_contains(0)[0]);
    load_ibf(index, "index_IBF_10");
    auto high = index.membership_agent();
    EXPECT_FALSE(high.bulk_contains(0)[0]);
    EXPECT_TRUE(high.bulk_contains(0x15'55'55'55'55'55'55'55ULL)[0]); // C repeated 31 times.
}

TEST_F(unitig_test, cli_reports_header_errors)
{
    std::ofstream{"invalid.fa"} << ">SRR1_0\n" << std::string(31, 'A') << '\n';
    for (std::string const mode : {"", "--ram"})
    {
        app_test_result result = execute_app("minimiser --unitigs -k 31 -w 31 -t 2", mode, "invalid.fa");
        EXPECT_FAILURE(result);
        EXPECT_EQ(result.err, "[Error] Missing unitig abundance (ka:f: or km:f:) in header: SRR1_0\n");
    }
}
