#include "utils.hpp"
#include <cassert>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

using namespace Eigen;


namespace utils 
{
    int factorial(int n)
    {
        int result = 1.;
        if (n != 0){
            for (int i= 1; i <= n; ++i){
                result *= i;
            }
        }
        return result;
    };

    void copyVec2Vec(const std::vector<double>& src, std::vector<double>& dst)
    {
        //assert (src.size() == dst.size());
        dst.resize(src.size());
        for (int i=0; i<src.size(); i++)
        {
            dst[i] = src[i];
        }
    };

    void copyEig2Vec(const VectorXd& src, std::vector<double>& dst)
    {
        assert (src.size() == dst.size());
        for (int i=0; i<src.size(); i++)
        {
            dst[i] = src(i);
        }
    };

    void copyEig2Arr(const VectorXd& src, double *dst)
    {
        for (int i=0; i<src.size(); i++){
            dst[i] = src(i);
        }
    };


    std::vector<double> copyEig2Vec(const VectorXd& src)
    {
        std::vector<double> dst(src.data(), src.data() + src.rows() * src.cols());
        return dst;
    };

   void copyMatvec2Matvec(const std::vector<std::vector<double>>& src, std::vector<std::vector<double>>& dst)
    {
        std::size_t d{src.size()};
        dst.resize(d);
        for (int i=0; i<d; i++)
        {
            dst[i].resize(src[i].size());
            copyVec2Vec(src[i], dst[i]);
        }
    };

    void writeVec2File(const std::vector<std::vector<double>> &data, const std::string &name)
    {
        std::ofstream file;
        std::string file_name;
        file_name = name + ".dat";
        file.open(file_name);
        std::size_t d{data.size()};
        std::size_t d2{data[0].size()};
        std::string line;
        if (file.is_open())
        {
            for (int i=0; i<d; i++)
            {
                line = "";
                for (int j=0; j<d2; j++)
                {
                   line += std::to_string(data[i][j]);
                   line += " "; 
                }
                line += "\n";
                file << line;
            }
            file.close();
        }
        else
        {
            std::cout << "unable to open file" << std::endl;
        }

    };

    void writeVec2File(const std::vector<std::vector<int>> &data, const std::string &name)
    {
        std::ofstream file;
        std::string file_name;
        file_name = name + ".dat";
        file.open(file_name);
        std::size_t d{data.size()};
        std::size_t d2{data[0].size()};
        std::string line;
        if (file.is_open())
        {
            for (int i=0; i<d; i++)
            {
                line = "";
                for (int j=0; j<d2; j++)
                {
                   line += std::to_string(data[i][j]);
                   line += " "; 
                }
                line += "\n";
                file << line;
            }
            file.close();
        }
        else
        {
            std::cout << "unable to open file" << std::endl;
        }

    };

    void writeMetadata2File(const std::string &problem_file, const std::string &name)
    {
        std::ifstream file;
        file.open(problem_file);
        std::ofstream metadata_file;
        std::string metadata_filename = name + ".meta";
        metadata_file.open(metadata_filename);
        std::string line;
        bool to_write{false};
        if (file.is_open() & metadata_file.is_open())
        {
            while (std::getline(file, line))
            {
                if (line.find("Problem") != std::string::npos)
                    to_write = true;
                if (to_write)
                {
                    if (line.find("int main") != std::string::npos)
                        to_write = false;
                    else
                       metadata_file << line << "\n";
                }
            }
            file.close();
            metadata_file.close();
        }
        else
        {
            std::cout << "unable to open file" << std::endl;
        }
    };
    
    std::string getDateString()
    {
        time_t now;
        char date_string[10];
        tm * now_tm;
        time(&now);
        now_tm = localtime(&now);
        strftime(date_string, 12, "%Y-%m-%d", now_tm);
        return date_string;
    };

    

}
