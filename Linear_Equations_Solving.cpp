/*--------------------------------------------------------------------------------------------------------------------------
* file name : Source.cpp
* Author : OUN
* Created on : July 30, 2019
* description: C++ program to solve linear equations using Gaussian elimination.
   The program reads first an integer n which is the number of equations.
   Modified to fix bugs and replace Cramer's rule with Gaussian elimination.
--------------------------------------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------INCLUDES------------------------------------------------------------*/
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>  // For fabs
#include <iomanip>  // For setprecision if needed
using namespace std;

/*-----------------------------------------------FUNCTIONS DEFINITIONS-----------------------------------------------------*/

// Improved maxmum function to find the highest variable index in an equation string
int maxmum(const string& xx) {
    int max_var = 0;
    size_t pos = 0;
    while ((pos = xx.find('x', pos)) != string::npos) {
        string no1 = xx.substr(pos + 1);
        int num1 = atoi(no1.c_str());  // Use atoi for integers
        if (num1 > max_var) {
            max_var = num1;
        }
        pos++;
    }
    return max_var;
}

// Gaussian elimination function
// Assumes square system (n equations, n variables)
// mat is [A | b] where A is n x n, b is n x 1 (augmented matrix)
// Solves Ax = b, stores solution in x[]
bool gaussian_elimination(float mat[100][101], int n, float x[100]) {
    // Forward elimination with partial pivoting
    for (int i = 0; i < n; ++i) {
        // Find pivot row
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(mat[k][i]) > fabs(mat[max_row][i])) {
                max_row = k;
            }
        }

        // Swap rows if necessary
        if (max_row != i) {
            for (int k = i; k <= n; ++k) {
                swap(mat[i][k], mat[max_row][k]);
            }
        }

        // Check for singular matrix
        if (fabs(mat[i][i]) < 1e-10) {
            return false;  // No unique solution
        }

        // Eliminate below pivot
        for (int k = i + 1; k < n; ++k) {
            float factor = mat[k][i] / mat[i][i];
            for (int j = i; j <= n; ++j) {
                mat[k][j] -= factor * mat[i][j];
            }
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        x[i] = mat[i][n];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= mat[i][j] * x[j];
        }
        x[i] /= mat[i][i];
    }

    return true;
}

int main() {
    int no_equ, x_pos = 0, z, max_var = 0;
    float num;
    string str, str2, no, no1;
    getline(cin, no);
    no_equ = atoi(no.c_str());  // Use atoi for integer

    vector<string> equ(no_equ);
    for (int i = 0; i < no_equ; i++) {
        getline(cin, equ[i]);
    }

    // First, find the global max variable index
    for (int i = 0; i < no_equ; i++) {
        max_var = max(max_var, maxmum(equ[i]));
    }

    // Now create the coefficient matrix ONCE
    vector<vector<float>> coeff(no_equ, vector<float>(max_var + 1, 0.0f));

    // Parse equations (this part is still fragile, but fixed some issues)
    for (int i = 0; i < no_equ; i++) {
        string eq = equ[i];  // Copy to modify if needed
        // Simple cleanup: remove spaces
        eq.erase(remove(eq.begin(), eq.end(), ' '), eq.end());

        size_t eq_pos = eq.find('=');
        if (eq_pos == string::npos) {
            // Invalid equation, but for now assume present
            continue;
        }

        // Parse left side
        string left = eq.substr(0, eq_pos);
        string right = eq.substr(eq_pos + 1);
        float constant = atof(right.c_str());
        coeff[i][0] = constant;  // Set constant directly

        // Parse terms in left
        size_t term_start = 0;
        while (term_start < left.length()) {
            size_t term_end = left.find_first_of("+-", term_start + 1);
            if (term_end == string::npos) term_end = left.length();

            string term = left.substr(term_start, term_end - term_start);
            float sign = 1.0f;
            if (term_start > 0 && left[term_start - 1] == '-') sign = -1.0f;

            size_t x_pos = term.find('x');
            float coef = 1.0f * sign;
            int var_idx = 0;
            if (x_pos != string::npos) {
                string coef_str = term.substr(0, x_pos);
                if (!coef_str.empty() && coef_str != "+" && coef_str != "-") {
                    coef = atof(coef_str.c_str()) * sign;
                }
                else if (coef_str == "-") {
                    coef = -1.0f;
                }
                string var_str = term.substr(x_pos + 1);
                var_idx = atoi(var_str.c_str());
            }
            else {
                // Constant on left, but since we move to right, subtract
                coef = atof(term.c_str()) * sign;
                coeff[i][0] -= coef;
                coef = 0;  // Not a variable
            }

            if (var_idx > 0) {
                coeff[i][var_idx] = coef;
            }

            term_start = term_end;
        }
    }

    //////////////////////////////////
    string operation;
    getline(cin, operation);

    /*----------------------------------------------- Level 1 -----------------------------------------------------*/
    if (operation == "num_vars") {
        cout << max_var;
    }
    else if (operation.substr(0, 8) == "equation") {
        string no_of_equ = operation.substr(8);
        int number = atoi(no_of_equ.c_str()) - 1;
        int max2 = maxmum(equ[number]);
        bool first = true;
        for (int i = 1; i <= max2; i++) {
            if (coeff[number][i] == 0) continue;
            if (!first && coeff[number][i] > 0) cout << "+";
            cout << coeff[number][i] << "x" << i;
            first = false;
        }
        cout << "=" << coeff[number][0];
    }
    else if (operation.substr(0, 8) == "column x") {
        string no_of_equ = operation.substr(8);
        int number = atoi(no_of_equ.c_str());
        for (int i = 0; i < no_equ; i++) {
            cout << coeff[i][number] << endl;
        }
    }
    /*----------------------------------------------- Level 2 -----------------------------------------------------*/
    else if (operation.substr(0, 3) == "add") {
        size_t space_pos = operation.find(' ', 4);
        string equ_1 = operation.substr(3, space_pos - 3);
        string equ_2 = operation.substr(space_pos + 1);
        int equ1 = atoi(equ_1.c_str());
        int equ2 = atoi(equ_2.c_str());
        int max1 = maxmum(equ[equ1 - 1]);
        int max2 = maxmum(equ[equ2 - 1]);
        int local_max = max(max1, max2);
        vector<float> result(local_max + 1, 0.0f);
        for (int i = 0; i <= local_max; i++) {
            result[i] = coeff[equ1 - 1][i] + coeff[equ2 - 1][i];
        }
        bool first = true;
        for (int i = 1; i <= local_max; i++) {
            if (result[i] == 0) continue;
            if (!first && result[i] > 0) cout << "+";
            cout << result[i] << "x" << i;
            first = false;
        }
        cout << "=" << result[0];
    }
    else if (operation.substr(0, 8) == "subtract") {
        size_t space_pos = operation.find(' ', 9);
        string equ_1 = operation.substr(8, space_pos - 8);
        string equ_2 = operation.substr(space_pos + 1);
        int equ1 = atoi(equ_1.c_str());
        int equ2 = atoi(equ_2.c_str());
        int max1 = maxmum(equ[equ1 - 1]);
        int max2 = maxmum(equ[equ2 - 1]);
        int local_max = max(max1, max2);
        vector<float> result(local_max + 1, 0.0f);
        for (int i = 0; i <= local_max; i++) {
            result[i] = coeff[equ1 - 1][i] - coeff[equ2 - 1][i];
        }
        bool first = true;
        for (int i = 1; i <= local_max; i++) {
            if (result[i] == 0) continue;
            if (!first && result[i] > 0) cout << "+";
            cout << result[i] << "x" << i;
            first = false;
        }
        cout << "=" << result[0];
    }
    else if (operation.substr(0, 12) == "substitute x") {
        size_t space1 = operation.find(' ', 13);
        size_t space2 = operation.find(' ', space1 + 1);
        string column_ = operation.substr(12, space1 - 12);
        string equ_1 = operation.substr(space1 + 1, space2 - space1 - 1);
        string equ_2 = operation.substr(space2 + 1);
        int column = atoi(column_.c_str());
        int equ1 = atoi(equ_1.c_str());
        int equ2 = atoi(equ_2.c_str());
        int local_max = max_var;  // Assume global max
        vector<float> new_equ(local_max + 1, 0.0f);
        vector<float> result(local_max + 1, 0.0f);
        float pivot = coeff[equ2 - 1][column];
        if (fabs(pivot) < 1e-10) {
            cout << "Division by zero in substitution";
            return 1;
        }
        float factor = -coeff[equ1 - 1][column] / pivot;
        for (int i = 0; i <= local_max; i++) {
            new_equ[i] = coeff[equ2 - 1][i] * factor;
        }
        for (int i = 0; i <= local_max; i++) {
            result[i] = new_equ[i] + coeff[equ1 - 1][i];
        }
        bool first = true;
        for (int i = 1; i <= local_max; i++) {
            if (result[i] == 0) continue;
            if (!first && result[i] > 0) cout << "+";
            cout << result[i] << "x" << i;
            first = false;
        }
        cout << "=" << result[0];
    }
    /*----------------------------------------------- Level 3 -----------------------------------------------------*/
    // Removed broken determinant functions
    // Replaced with Gaussian for "solve"
    else if (operation == "solve") {
        if (no_equ != max_var) {
            cout << "System is not square (equations != variables). Cannot solve uniquely.";
            return 1;
        }
        int n = no_equ;  // Square system
        float mat[100][101];  // Augmented [A | b]
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                mat[i][j] = coeff[i][j + 1];  // A
            }
            mat[i][n] = coeff[i][0];  // b
        }

        float solution[100];
        if (gaussian_elimination(mat, n, solution)) {
            for (int k = 0; k < n; k++) {
                cout << "x" << (k + 1) << "=" << solution[k] << endl;
            }
        }
        else {
            cout << "No unique solution (singular matrix or inconsistent system).";
        }
    }


    return 0;
}