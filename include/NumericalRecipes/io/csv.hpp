
#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <NumericalRecipes/types.hpp>


namespace NumericalRecipes {
namespace io {

    namespace {
        using namespace types;
    }

    class CSVReaderBase{
        public:
            std::string file_name;
            std::string delimiter;
            std::string current_line;

            int n_skip_lines;
        public:
            CSVReaderBase(
                std::string file_name_,
                std::string delimiter_=",",
                int n_skip_lines_=0)
            :
                file_name(file_name_),
                delimiter(delimiter_),
                n_skip_lines(n_skip_lines_),
                current_line("")
            {}
    };

    class TableReader: public CSVReaderBase{

        public:

            public:
                typedef typename CSVReaderBase Super;

            public:

                using Super::CSVReaderBase;

            unordered_map<string, VectorDouble> get_data(){
                std::ifstream file(file_name);

                getline(file, current_line);
                boost::algorithm::to_lower(current_line);
                std::vector<std::string> category_vector;
                boost::algorithm::split(category_vector, current_line, boost::is_any_of(delimiter));

                unordered_map<string,VectorDouble> data_umap;

                for (int i=0; i<category_vector.size(); i++)
                    data_umap[category_vector[i]] = VectorDouble();

                while (getline(file, current_line)){
                    std::vector<std::string> vector_string;
                    boost::algorithm::split(vector_string, current_line, boost::is_any_of(delimiter));
                    for (int i=0; i<category_vector.size(); i++)
                        data_umap[category_vector[i]].push_back(std::stod(vector_string[i]));
               }

                return data_umap;
            }
    };
}
}