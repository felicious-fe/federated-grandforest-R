/*------------------------------------------------------------------------------
 This file is part of Grand Forest.

 Grand Forest incorporates work that is part of Ranger.
 See the original license notice below.
------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck

 http://www.imbs-luebeck.de
 #-------------------------------------------------------------------------------*/

#include <RcppEigen.h>
#include <vector>
#include <sstream>

#include "globals.h"
#include "Forest.h"
#include "ForestClassification.h"
#include "ForestRegression.h"
#include "ForestSurvival.h"
#include "ForestProbability.h"
#include "Data.h"
#include "DataChar.h"
#include "DataDouble.h"
#include "DataFloat.h"
#include "DataSparse.h"

class List;

class List;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List grandforestCpp(
        uint treetype, std::string dependent_variable_name,
        Rcpp::NumericMatrix input_data, Rcpp::NumericMatrix graph_data,
        std::vector <std::string> variable_names, uint mtry, uint num_trees, bool verbose,
        uint seed, uint num_threads, bool write_forest, uint importance_mode_r, uint subgraph_mode_r,
        uint min_node_size, std::vector <std::vector<double>> &split_select_weights, bool use_split_select_weights,
        std::vector <std::string> &always_split_variable_names, bool use_always_split_variable_names,
        std::string status_variable_name, bool prediction_mode, Rcpp::List loaded_forest, Rcpp::RawMatrix snp_data,
        bool sample_with_replacement, bool probability, std::vector <std::string> &unordered_variable_names,
        bool use_unordered_variable_names, bool save_memory, uint splitrule_r,
        std::vector<double> &case_weights, bool use_case_weights, bool predict_all,
        bool keep_inbag, double sample_fraction, double alpha, double minprop, bool holdout, uint prediction_type_r,
        uint num_random_splits, bool random_root, Eigen::SparseMatrix<double> sparse_data, bool use_sparse_data
) {
    Rcpp::List result;
    Forest *forest = 0;
    Data *data = 0;
    Graph *graph = 0;
    try {
        // Empty split select weights and always split variables if not used
        if (!use_split_select_weights) {
            split_select_weights.clear();
        }
        if (!use_always_split_variable_names) {
            always_split_variable_names.clear();
        }
        if (!use_unordered_variable_names) {
            unordered_variable_names.clear();
        }
        if (!use_case_weights) {
            case_weights.clear();
        }

        std::ostream *verbose_out;
        if (verbose) {
            verbose_out = &Rcpp::Rcout;
        } else {
            verbose_out = new std::stringstream;
        }

        size_t num_rows;
        size_t num_cols;
        if (use_sparse_data) {
            num_rows = sparse_data.rows();
            num_cols = sparse_data.cols();
        } else {
            num_rows = input_data.nrow();
            num_cols = input_data.ncol();
        }

        // Initialize data
        if (use_sparse_data) {
            data = new DataSparse(&sparse_data, variable_names, num_rows, num_cols);
        } else {
            data = new DataDouble(input_data.begin(), variable_names, num_rows, num_cols);
        }

        // If there is snp data, add it
        if (snp_data.nrow() > 1) {
            data->addSnpData(snp_data.begin(), snp_data.ncol());
        }

        // Initialize graph data
        if (graph_data.nrow() > 1) {
            std::vector <std::pair<size_t, size_t>> edges;
            for (size_t i = 0; i < graph_data.nrow(); ++i) {
                edges.emplace_back(graph_data(i, 0), graph_data(i, 1));
            }
            graph = new Graph(data->getNumCols(), edges);
        }

        switch (treetype) {
            case TREE_CLASSIFICATION:
                if (probability) {
                    forest = new ForestProbability;
                } else {
                    forest = new ForestClassification;
                }
                break;
            case TREE_REGRESSION:
                forest = new ForestRegression;
                break;
            case TREE_SURVIVAL:
                forest = new ForestSurvival;
                break;
            case TREE_PROBABILITY:
                forest = new ForestProbability;
                break;
        }

        ImportanceMode importance_mode = (ImportanceMode) importance_mode_r;
        SubgraphMode subgraph_mode = (SubgraphMode) subgraph_mode_r;
        SplitRule splitrule = (SplitRule) splitrule_r;
        PredictionType prediction_type = (PredictionType) prediction_type_r;

        // Init Grand Forest
        forest->initR(dependent_variable_name, data, graph, mtry, num_trees, verbose_out, seed, num_threads,
                      importance_mode, subgraph_mode, min_node_size, split_select_weights, always_split_variable_names,
                      status_variable_name, prediction_mode, sample_with_replacement, unordered_variable_names,
                      save_memory, splitrule, case_weights, predict_all, keep_inbag, sample_fraction, alpha, minprop,
                      holdout, prediction_type, num_random_splits, random_root);

        // Load forest object if in prediction mode
        if (prediction_mode) {
            size_t dependent_varID = loaded_forest["dependent.varID"];
            //size_t num_trees = loaded_forest["num.trees"];
            std::vector < std::vector < std::vector < size_t>> > child_nodeIDs = loaded_forest["child.nodeIDs"];
            std::vector <std::vector<size_t>> split_varIDs = loaded_forest["split.varIDs"];
            std::vector <std::vector<double>> split_values = loaded_forest["split.values"];
            std::vector<bool> is_ordered = loaded_forest["is.ordered"];

            if (treetype == TREE_CLASSIFICATION) {
                std::vector<double> class_values = loaded_forest["class.values"];
                ((ForestClassification *) forest)->loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs,
                                                              split_values, class_values, is_ordered);
            } else if (treetype == TREE_REGRESSION) {
                ((ForestRegression *) forest)->loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs,
                                                          split_values,
                                                          is_ordered);
            } else if (treetype == TREE_SURVIVAL) {
                size_t status_varID = loaded_forest["status.varID"];
                std::vector < std::vector < std::vector < double>> > chf = loaded_forest["chf"];
                std::vector<double> unique_timepoints = loaded_forest["unique.death.times"];
                ((ForestSurvival *) forest)->loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs,
                                                        split_values,
                                                        status_varID, chf, unique_timepoints, is_ordered);
            } else if (treetype == TREE_PROBABILITY) {
                std::vector<double> class_values = loaded_forest["class.values"];
                std::vector < std::vector < std::vector < double>>>terminal_class_counts =
                                                                           loaded_forest["terminal.class.counts"];
                ((ForestProbability *) forest)->loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs,
                                                           split_values,
                                                           class_values, terminal_class_counts, is_ordered);
            }
        }

        // Run Grand Forest
        forest->run(false);

        if (use_split_select_weights && importance_mode != IMP_NONE) {
            *verbose_out
                    << "Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights."
                    << std::endl;
        }

        // Use first non-empty dimension of predictions
        const std::vector <std::vector<std::vector < double>>>&predictions = forest->getPredictions();
        if (predictions.size() == 1) {
            if (predictions[0].size() == 1) {
                result.push_back(forest->getPredictions()[0][0], "predictions");
            } else {
                result.push_back(forest->getPredictions()[0], "predictions");
            }
        } else {
            result.push_back(forest->getPredictions(), "predictions");
        }

        // Return output
        result.push_back(forest->getNumTrees(), "num.trees");
        result.push_back(forest->getNumIndependentVariables(), "num.independent.variables");
        if (treetype == TREE_SURVIVAL) {
            ForestSurvival *temp = (ForestSurvival *) forest;
            result.push_back(temp->getUniqueTimepoints(), "unique.death.times");
        }
        if (!verbose) {
            std::stringstream temp;
            temp << verbose_out->rdbuf();
            result.push_back(temp.str(), "log");
        }
        if (!prediction_mode) {
            result.push_back(forest->getMtry(), "mtry");
            result.push_back(forest->getMinNodeSize(), "min.node.size");
            if (importance_mode != IMP_NONE) {
                result.push_back(forest->getVariableImportance(), "variable.importance");
            }
            result.push_back(forest->getOverallPredictionError(), "prediction.error");
            result.push_back(forest->getVariableFrequency(), "variable.frequency");
        }

        if (keep_inbag) {
            result.push_back(forest->getInbagCounts(), "inbag.counts");
        }

        // Save forest if needed
        if (write_forest) {
            Rcpp::List forest_object;
            forest_object.push_back(forest->getDependentVarId(), "dependent.varID");
            forest_object.push_back(forest->getNumTrees(), "num.trees");
            forest_object.push_back(forest->getChildNodeIDs(), "child.nodeIDs");
            forest_object.push_back(forest->getSplitVarIDs(), "split.varIDs");
            forest_object.push_back(forest->getSplitValues(), "split.values");
            forest_object.push_back(forest->getIsOrderedVariable(), "is.ordered");

            if (treetype == TREE_CLASSIFICATION) {
                ForestClassification *temp = (ForestClassification *) forest;
                forest_object.push_back(temp->getClassValues(), "class.values");
            } else if (treetype == TREE_PROBABILITY) {
                ForestProbability *temp = (ForestProbability *) forest;
                forest_object.push_back(temp->getClassValues(), "class.values");
                forest_object.push_back(temp->getTerminalClassCounts(), "terminal.class.counts");
            } else if (treetype == TREE_SURVIVAL) {
                ForestSurvival *temp = (ForestSurvival *) forest;
                forest_object.push_back(temp->getStatusVarId(), "status.varID");
                forest_object.push_back(temp->getChf(), "chf");
                forest_object.push_back(temp->getUniqueTimepoints(), "unique.death.times");
            }
            result.push_back(forest_object, "forest");
        }

        delete forest;
        delete data;
        if (graph != 0) delete graph;
    } catch (std::exception &e) {
        if (strcmp(e.what(), "User interrupt.") != 0) {
            Rcpp::Rcerr << "Error: " << e.what() << " Grand Forest will EXIT now." << std::endl;
        }
        delete forest;
        delete data;
        if (graph != 0) delete graph;
        return result;
    }

    return result;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List grandforest_sum_modelsCpp(Rcpp::List rmodel1, Rcpp::List rmodel2,
                                     uint treetype, std::string dependent_variable_name, Rcpp::NumericMatrix input_data1,
                                     Rcpp::NumericMatrix input_data2,std::vector <std::string> variable_names1,
                                     std::vector <std::string> variable_names2,
                                     Rcpp::NumericMatrix graph_data, uint mtry, size_t num_trees1, size_t num_trees2, bool verbose, uint seed,
                                     uint num_threads, uint importance_mode_r, uint subgraph_mode_r, uint min_node_size,
                                     std::vector <std::vector<double>> &split_select_weights,
                                     bool use_split_select_weights,
                                     std::vector <std::string> &always_split_variable_names,
                                     bool use_always_split_variable_names, std::string status_variable_name,
                                     bool prediction_mode, bool sample_with_replacement, bool probability,
                                     std::vector <std::string> &unordered_variable_names,
                                     bool use_unordered_variable_names, bool save_memory, uint splitrule_r,
                                     std::vector<double> &case_weights, bool use_case_weights, bool predict_all,
                                     bool keep_inbag, double sample_fraction, double alpha, double minprop,
                                     bool holdout, uint prediction_type_r, uint num_random_splits, bool random_root)
{
    Rcpp::Rcout << "Summing models..." << std::endl;

    Rcpp::List result;
    Forest *forest1 = nullptr;
    Forest *forest2 = nullptr;
    Data *data1 = nullptr;
    Data *data2 = nullptr;
    Graph *graph = nullptr;

    try {
        Rcpp::List rforest1 = rmodel1["forest"];
        Rcpp::List rforest2 = rmodel2["forest"];

        // Empty split select weights and always split variables if not used
        if (!use_split_select_weights) {
            split_select_weights.clear();
        }
        if (!use_always_split_variable_names) {
            always_split_variable_names.clear();
        }
        if (!use_unordered_variable_names) {
            unordered_variable_names.clear();
        }
        if (!use_case_weights) {
            case_weights.clear();
        }

        auto importance_mode = (ImportanceMode) importance_mode_r;
        auto subgraph_mode = (SubgraphMode) subgraph_mode_r;
        auto splitrule = (SplitRule) splitrule_r;
        auto prediction_type = (PredictionType) prediction_type_r;


        std::ostream *verbose_out;
        if (verbose) {
            verbose_out = &Rcpp::Rcout;
        } else {
            verbose_out = new std::stringstream;
        }

        // create dummy Data objects to correctly initialize the forest objects
        size_t num_rows1;
        size_t num_cols1;
        num_rows1 = input_data1.nrow();
        num_cols1 = input_data1.ncol();
        size_t num_rows2;
        size_t num_cols2;
        num_rows2 = input_data2.nrow();
        num_cols2 = input_data2.ncol();
        data1 = new DataDouble(input_data1.begin(), variable_names1, num_rows1, num_cols1);
        data2 = new DataDouble(input_data2.begin(), variable_names2, num_rows2, num_cols2);

        switch (treetype) {
            case TREE_CLASSIFICATION:
                if (probability) {
                    forest1 = new ForestProbability;
                    forest2 = new ForestProbability;
                } else {
                    forest1 = new ForestClassification;
                    forest2 = new ForestClassification;
                }
                break;
            case TREE_REGRESSION:
                forest1 = new ForestRegression;
                forest2 = new ForestRegression;
                break;
            case TREE_SURVIVAL:
                forest1 = new ForestSurvival;
                forest2 = new ForestSurvival;
                break;
            case TREE_PROBABILITY:
                forest1 = new ForestProbability;
                forest2 = new ForestProbability;
                break;
            default:
                Rcpp::Rcout << "WRONG TREETYPE" << std::endl;
                break;
        }

        forest1->initR(dependent_variable_name, data1, graph, mtry, num_trees1, verbose_out, seed, num_threads,
                       importance_mode, subgraph_mode, min_node_size, split_select_weights, always_split_variable_names,
                       status_variable_name, prediction_mode, sample_with_replacement, unordered_variable_names,
                       save_memory, splitrule, case_weights, predict_all, keep_inbag, sample_fraction, alpha, minprop,
                       holdout, prediction_type, num_random_splits, random_root);

        forest2->initR(dependent_variable_name, data2, graph, mtry, num_trees2, verbose_out, seed, num_threads,
                      importance_mode, subgraph_mode, min_node_size, split_select_weights, always_split_variable_names,
                      status_variable_name, prediction_mode, sample_with_replacement, unordered_variable_names,
                      save_memory, splitrule, case_weights, predict_all, keep_inbag, sample_fraction, alpha, minprop,
                      holdout, prediction_type, num_random_splits, random_root);


        // load forest1
        size_t dependent_varID1 = rforest1["dependent.varID"];
        std::vector < std::vector < std::vector < size_t>> > child_nodeIDs1 = rforest1["child.nodeIDs"];
        std::vector <std::vector<size_t>> split_varIDs1 = rforest1["split.varIDs"];
        std::vector <std::vector<double>> split_values1 = rforest1["split.values"];
        std::vector<bool> is_ordered1 = rforest1["is.ordered"];
        std::vector<double> variable_importance1 = rmodel1["variable.importance"];
        std::vector<std::string> independent_variable_names1 = rforest1["independent.variable.names"];

        std::map<std::string, double> variable_importance_map1;
        assert(independent_variable_names1.size() == variable_importance1.size());
        for (size_t i = 0; i < independent_variable_names1.size(); ++i)
            variable_importance_map1[independent_variable_names1[i]] = variable_importance1[i];

        // load forest2
        size_t dependent_varID2 = rforest2["dependent.varID"];
        std::vector < std::vector < std::vector < size_t>> > child_nodeIDs2 = rforest2["child.nodeIDs"];
        std::vector <std::vector<size_t>> split_varIDs2 = rforest2["split.varIDs"];
        std::vector <std::vector<double>> split_values2 = rforest2["split.values"];
        std::vector<bool> is_ordered2 = rforest2["is.ordered"];
        std::vector<double> variable_importance2 = rmodel2["variable.importance"];
        std::vector<std::string> independent_variable_names2 = rforest2["independent.variable.names"];

        std::map<std::string, double> variable_importance_map2;
        assert(independent_variable_names2.size() == variable_importance2.size());
        for (size_t i = 0; i < independent_variable_names2.size(); ++i)
            variable_importance_map2[independent_variable_names2[i]] = variable_importance2[i];


        if (treetype == TREE_CLASSIFICATION) {
            if (probability) {
                std::vector<double> class_values1 = rforest1["class.values"];
                ((ForestClassification *) forest1)->loadForest(dependent_varID1, num_trees1, child_nodeIDs1,
                                                               split_varIDs1,
                                                               split_values1, class_values1, is_ordered1);
                std::vector<double> class_values2 = rforest2["class.values"];
                ((ForestClassification *) forest2)->loadForest(dependent_varID2, num_trees2, child_nodeIDs2,
                                                               split_varIDs2,
                                                               split_values2, class_values2, is_ordered2);
            } else {
                ((ForestRegression *) forest1)->loadForest(dependent_varID1, num_trees1, child_nodeIDs1, split_varIDs1,
                                                           split_values1, is_ordered1);
                ((ForestRegression *) forest2)->loadForest(dependent_varID2, num_trees2, child_nodeIDs2, split_varIDs2,
                                                           split_values2, is_ordered2);
            }
        } else if (treetype == TREE_REGRESSION) {
            ((ForestRegression *) forest1)->loadForest(dependent_varID1, num_trees1, child_nodeIDs1, split_varIDs1, split_values1, is_ordered1);
            ((ForestRegression *) forest2)->loadForest(dependent_varID2, num_trees2, child_nodeIDs2, split_varIDs2, split_values2, is_ordered2);
        } else if (treetype == TREE_SURVIVAL) {
            size_t status_varID1 = rforest1["status.varID"];
            std::vector < std::vector < std::vector < double>> > chf1 = rforest1["chf"];
            std::vector<double> unique_timepoints1 = rforest1["unique.death.times"];
            ((ForestSurvival *) forest1)->loadForest(dependent_varID1, num_trees1, child_nodeIDs1, split_varIDs1,
                                                     split_values1,
                                                     status_varID1, chf1, unique_timepoints1, is_ordered1);
            size_t status_varID2 = rforest2["status.varID"];
            std::vector < std::vector < std::vector < double>> > chf2 = rforest2["chf"];
            std::vector<double> unique_timepoints2 = rforest2["unique.death.times"];
            ((ForestSurvival *) forest2)->loadForest(dependent_varID2, num_trees2, child_nodeIDs2, split_varIDs2,
                                                     split_values2,
                                                     status_varID2, chf2, unique_timepoints2, is_ordered2);
        } else if (treetype == TREE_PROBABILITY) {
            std::vector<double> class_values1 = rforest1["class.values"];
            std::vector < std::vector < std::vector <
            double>>>terminal_class_counts1 = rforest1["terminal.class.counts"];
            ((ForestProbability *) forest1)->loadForest(dependent_varID1, num_trees1, child_nodeIDs1, split_varIDs1,
                                                        split_values1,
                                                        class_values1, terminal_class_counts1, is_ordered1);
            std::vector<double> class_values2 = rforest2["class.values"];
            std::vector < std::vector < std::vector <
            double>>>terminal_class_counts2 = rforest2["terminal.class.counts"];
            ((ForestProbability *) forest2)->loadForest(dependent_varID2, num_trees2, child_nodeIDs2, split_varIDs2,
                                                        split_values2,
                                                        class_values2, terminal_class_counts2, is_ordered2);
        }

        Rcpp::Rcout << "CREATED FOREST OBJECTS" << std::endl;

        // Sum forests
        std::vector<Tree *> forest2_trees = forest2->getTrees();
        forest1->addTrees(forest2_trees);
        forest1->computeVariableFrequency();

        Rcpp::Rcout << "COMBINED TREE SETS" << std::endl;

        // combine variable importances
        std::map<std::string, double> combined_variable_importance_map;
        std::vector<std::string> combined_variable_names;
        std::set_intersection(independent_variable_names1.begin(), independent_variable_names1.end(),
                              independent_variable_names2.begin(), independent_variable_names2.end(),
                              std::back_inserter(combined_variable_names));

        std::vector<int>::iterator it, ls;

        for (const auto name : combined_variable_names) {
            Rcpp::Rcout << name << ";";
        }
        Rcpp::Rcout << std::endl;

        size_t max_importance = 0;
        for (const auto& name : combined_variable_names) {
            auto pos1 = variable_importance_map1.find(name);
            double value1;
            if (pos1 == variable_importance_map1.end()) {
                value1 = 0;
            } else {
                value1 = pos1->second;
            }
            auto pos2 = variable_importance_map2.find(name);
            double value2;
            if (pos2 == variable_importance_map2.end()) {
                value2 = 0;
            } else {
                value2 = pos2->second;
            }

            size_t importance = value1 * (double)rmodel1["num.samples"] + value2 * (double)rmodel2["num.samples"];

            if (importance > max_importance) {
                max_importance = importance;
            }

            combined_variable_importance_map[name] = importance;
        }

        Rcpp::DoubleVector combined_variable_importance = {};
        Rcpp::StringVector combined_variable_names_rcpp = {};
        for (const auto& name : combined_variable_names) {
            double importance = combined_variable_importance_map.find(name)->second;
            combined_variable_importance.push_back(importance / max_importance);
            Rcpp::String str = {name};
            combined_variable_names_rcpp.push_back(str);
            //combined_variable_importance_map[name] = importance / max_importance;
        }

        combined_variable_importance.names() = combined_variable_names_rcpp;

        Rcpp::Rcout << "SUMMED FOREST OBJECTS" << std::endl;

        // Return output
        result.push_back(forest1->getNumTrees(), "num.trees");
        result.push_back(forest1->getNumIndependentVariables(), "num.independent.variables");

        if (treetype == TREE_SURVIVAL) {
            ForestSurvival *temp = (ForestSurvival *) forest1;
            result.push_back(temp->getUniqueTimepoints(), "unique.death.times");
        }
        if (!verbose) {
            std::stringstream temp;
            temp << verbose_out->rdbuf();
            result.push_back(temp.str(), "log");
        }
        result.push_back(forest1->getMtry(), "mtry");
        result.push_back(forest1->getMinNodeSize(), "min.node.size");
        result.push_back(combined_variable_importance, "variable.importance");
        result.push_back(forest1->getOverallPredictionError(), "prediction.error");
        result.push_back(forest1->getVariableFrequency(), "variable.frequency");

        if (keep_inbag) {
            result.push_back(forest1->getInbagCounts(), "inbag.counts");
        }

        // Save forest
        Rcpp::List forest_object;
        forest_object.push_back(forest1->getDependentVarId(), "dependent.varID");
        forest_object.push_back(forest1->getNumTrees(), "num.trees");
        forest_object.push_back(forest1->getChildNodeIDs(), "child.nodeIDs");
        forest_object.push_back(forest1->getSplitVarIDs(), "split.varIDs");
        forest_object.push_back(forest1->getSplitValues(), "split.values");
        forest_object.push_back(forest1->getIsOrderedVariable(), "is.ordered");

        if (treetype == TREE_CLASSIFICATION) {
            ForestClassification *temp = (ForestClassification *) forest1;
            forest_object.push_back(temp->getClassValues(), "class.values");
        } else if (treetype == TREE_PROBABILITY) {
            ForestProbability *temp = (ForestProbability *) forest1;
            forest_object.push_back(temp->getClassValues(), "class.values");
            forest_object.push_back(temp->getTerminalClassCounts(), "terminal.class.counts");
        } else if (treetype == TREE_SURVIVAL) {
            ForestSurvival *temp = (ForestSurvival *) forest1;
            forest_object.push_back(temp->getStatusVarId(), "status.varID");
            forest_object.push_back(temp->getChf(), "chf");
            forest_object.push_back(temp->getUniqueTimepoints(), "unique.death.times");
        }
        result.push_back(forest_object, "forest");

        Rcpp::Rcout << "CLEANING UP..." << std::endl;

        delete forest1;
        // Forest 2 does not need to be deleted, the destructor of Forest 1 already deletes the Trees
        //delete forest2;
        delete data1;
        delete data2;
    } catch (std::exception &e) {
        if (strcmp(e.what(), "User interrupt.") != 0) {
            Rcpp::Rcerr << "Error: " << e.what() << " Grand Forest will EXIT now." << std::endl;
        }
        delete forest1;
        // Forest 2 does not need to be deleted, the destructor of Forest 1 already deletes the Trees
        //delete forest2;
        return result;
    }

    return result;

}
