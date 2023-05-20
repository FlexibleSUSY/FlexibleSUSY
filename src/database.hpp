// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include <string>
#include <vector>
#include <Eigen/Core>

struct sqlite3;

namespace flexiblesusy {
namespace database {

class Database {
public:
   Database(const std::string& file_name);
   Database(const Database&) = delete;
   Database(Database&&) = delete;
   ~Database();

   /// insert a row of doubles into a table
   void insert(const std::string&, const std::vector<std::string>&, const Eigen::ArrayXd&);

   /// extract a row of doubles from a table
   Eigen::ArrayXd extract(const std::string&, long long);

private:
   using TCallback = int (*)(void*, int, char**, char**);

   sqlite3* db{nullptr}; ///< pointer to database object

   void execute(const std::string&);
   void execute(const std::string&, TCallback, void*);
   template <typename T>
   void create_table(const std::string&, const std::vector<std::string>&);
};

} // namespace database
} // namespace flexiblesusy
