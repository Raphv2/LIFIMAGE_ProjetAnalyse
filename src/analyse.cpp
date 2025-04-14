/***************************************************************************************************
 * @file  analyse.cpp
 * @brief Contains the main program of the project
 **************************************************************************************************/
#include <opencv4/opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <cstdlib> 
#include <ctime>   
#include <cmath>

/**
 * @brief Calcule la couleur dominante dans une image HSV.
 * 
 * @param imageHSV Image convertie en espace HSV.
 * @return cv::Scalar La couleur dominante (H, S, V).
 */
cv::Scalar getDominantColor(cv::Mat& imageHSV) {
    std::vector<cv::Mat> hsvChannels;
    cv::split(imageHSV, hsvChannels);

    int histSize = 180; 
    float range[] = {0, 180}; 
    const float* histRange = {range};
    cv::Mat hist;

    cv::calcHist(&hsvChannels[0], 1, 0, cv::Mat(), hist, 1, &histSize, &histRange, true, false);

    double maxVal = 0;
    cv::Point maxLoc; 
    cv::minMaxLoc(hist, 0, &maxVal, 0, &maxLoc);
    
    int dominantHue = maxLoc.y; 

    double avgSaturation = cv::mean(hsvChannels[1])[0];
    double avgValue = cv::mean(hsvChannels[2])[0];

    return cv::Scalar(dominantHue, avgSaturation, avgValue);  
}

/**
 * @brief Génère une couleur RGB aléatoire.
 * 
 * @return cv::Scalar Couleur aléatoire (B, G, R).
 */
cv::Scalar getRandomColor() {
return cv::Scalar(std::rand() % 256, std::rand() % 256, std::rand() % 256);
}

/**
 * @brief Simplifie les contours à l’aide de l’algorithme Douglas-Peucker.
 * 
 * @param contours Liste des contours d’entrée.
 * @param coef Coefficient déterminant la précision de la simplification.
 * @return std::vector<std::vector<cv::Point>> Contours simplifiés.
 */
std::vector<std::vector<cv::Point>> simplifiedContour( std::vector<std::vector<cv::Point>> contours, double coef){
    std::vector<std::vector<cv::Point>> result;

    for (size_t i = 0; i < contours.size(); ++i) {
        std::vector<cv::Point> simplified;
        
        double epsilon = coef * cv::arcLength(contours[i], true); 
        cv::approxPolyDP(contours[i], simplified, epsilon, true); 

        result.push_back(simplified);
    }
    return result;
}

/**
 * @brief Extrait les objets dans des vignettes à partir des contours et les enregistre.
 * 
 * @param contours Contours détectés.
 * @param contourImage Image contenant les objets.
 */
void vignetHandler(std::vector<std::vector<cv::Point>> contours, cv::Mat& contourImage){
    int objIndex = 0;
    for (const auto& contour : contours) {

        cv::Rect boundingBox = cv::boundingRect(contour);

        cv::Mat roi = contourImage(boundingBox);

        std::string filename = "data/object_" + std::to_string(objIndex) + ".png";
        cv::imwrite(filename, roi);
        objIndex++;
    }
}

/**
 * @brief Calcule l’angle entre trois points.
 * 
 * @param A Premier point.
 * @param B Point central.
 * @param C Troisième point.
 * @return double Angle en degrés.
 */
double Angle(cv::Point A, cv::Point B, cv::Point C) {
    cv::Point2f AB = B - A;
    cv::Point2f BC = C - B;

    double dot = AB.x * BC.x + AB.y * BC.y;
    double det = AB.x * BC.y - AB.y * BC.x;  

    double angleRad = std::atan2(det, dot); 
    double angleDeg = angleRad * 180.0 / CV_PI;

    return angleDeg; 
}

/**
 * @brief Dessine les lignes entre les points des contours et affiche les angles convexes.
 * 
 * @param contours Contours détectés.
 * @param contourImage Image où les lignes sont dessinées.
 */
void drawLineHandler(std::vector<std::vector<cv::Point>> contours, cv::Mat& contourImage){
    for (size_t i = 0; i < contours.size(); ++i) {
        double area = cv::contourArea(contours[i]);
        if (area < 1000) continue; 

        std::vector<cv::Point> approx;
        cv::approxPolyDP(contours[i], approx, 16.0, true);

        for (size_t j = 0; j < approx.size(); ++j) {

            cv::Point A = approx[(j - 1 + approx.size()) % approx.size()];
            cv::Point B = approx[j];
            cv::Point C = approx[(j + 1) % approx.size()];


            cv::line(contourImage, A, B, cv::Scalar(0, 0, 0), 2);

            double ang_contour = Angle(A, B, C);

            if (ang_contour > 30.0 ) {
                std::cout << "Concave à " << j << " (angle = " << ang_contour << "°)" << std::endl;
                cv::line(contourImage, B, B, cv::Scalar(0, 0, 255), 10); 
                cv::line(contourImage, A, C, cv::Scalar(255, 0, 0), 2);
            } 
            if (ang_contour < -30.0 ) {
                std::cout << "Convexe à " << j << " (angle = " << ang_contour << "°)" << std::endl;
                cv::line(contourImage, B, B, cv::Scalar(0, 255, 0), 10); 
            } 
        }
    }
}

/**
 * @brief Point d’entrée du programme.
 * 
 * Charge une image, détecte sa couleur dominante, crée un masque, détecte les contours,
 * dessine les contours simplifiés avec angles, extrait les vignettes et affiche les résultats.
 * 
 * @return int Code de sortie.
 */
int main() {
    cv::Mat inputImage = cv::imread("data/puzzle.jpg");
    if (inputImage.empty()) {
        std::cerr << "Erreur lors du chargement de l'image." << std::endl;
        return -1;
    }

    int width = inputImage.cols;
    int height = inputImage.rows;

    cv::Rect resize(10 , 10 , width * 0.9, height - 20); 

    cv::Mat croppedImage = inputImage(resize);



    cv::Mat imageHSV;
    cv::cvtColor(croppedImage, imageHSV, cv::COLOR_BGR2HSV);

    cv::Scalar dominantColor = getDominantColor(imageHSV);
    std::cout << "Teinte dominante : " << dominantColor[0] << std::endl;

    float hueRange = 6.5; 

    cv::Scalar lowerBound(dominantColor[0] - hueRange, dominantColor[1], dominantColor[2]);
    cv::Scalar upperBound(dominantColor[0] + hueRange, 255, 255);

    cv::Mat mask;
    cv::inRange(imageHSV, lowerBound, upperBound, mask);
    cv::bitwise_not(mask, mask);

    cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(5, 5));
    cv::erode(mask, mask, kernel);
    cv::dilate(mask, mask, kernel);

    std::vector<std::vector<cv::Point>> contours;
    std::vector<cv::Vec4i> hierarchy;
    cv::findContours(mask, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    cv::Mat contourImage = cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

    for (size_t i = 0; i < contours.size(); ++i) {

        cv::Mat pieceMask = cv::Mat::zeros(croppedImage.size(), CV_8UC1);
        cv::drawContours(pieceMask, contours, static_cast<int>(i), cv::Scalar(255), cv::FILLED);

        cv::Rect boundingBox = cv::boundingRect(contours[i]);

        cv::Mat piece;
        croppedImage.copyTo(piece, pieceMask); 

        cv::Mat croppedPiece = piece(boundingBox);
        cv::Mat croppedMask = pieceMask(boundingBox);

        croppedPiece.copyTo(contourImage(boundingBox), croppedMask);
    }

    drawLineHandler(contours, contourImage);

    vignetHandler(contours, contourImage);

    cv::imshow("Image Originale", croppedImage);
    cv::imshow("Masque basé sur la couleur dominante (orange)", mask);
    cv::imshow("Contours trouvés", contourImage);

    cv::waitKey(0);

    return 0;
 }
 

 