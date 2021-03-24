compilation:

from test folder

g++ -I ../src kde_test.cpp

https://stackoverflow.com/questions/25274312/is-it-a-good-practice-to-define-c-functions-inside-header-files


g++ -I ../src kde_test.cpp ../src/kde.cpp ../src/utils.cpp -o kde_test


include nlopt
g++ -I ../src optimizer_test.cpp -o optimizer_test -lnlopt
