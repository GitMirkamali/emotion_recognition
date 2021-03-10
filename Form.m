function varargout = Form(varargin)
% FORM M-file for Form.fig
%      FORM, by itself, creates a new FORM or raises the existing
%      singleton*.
%
%      H = FORM returns the handle to a new FORM or the handle to
%      the existing singleton*.
%
%      FORM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FORM.M with the given input arguments.
%
%      FORM('Property','Value',...) creates a new FORM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Form_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Form_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Form

% Last Modified by GUIDE v2.5 11-Aug-2010 20:08:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Form_OpeningFcn, ...
                   'gui_OutputFcn',  @Form_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Form is made visible.
function Form_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Form (see VARARGIN)

% Choose default command line output for Form
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Form wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Form_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in bttOpenImage.
function bttOpenImage_Callback(hObject, eventdata, handles)
% hObject    handle to bttOpenImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Face;

[FileName, PathName] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'},'Select an Image');
fpath = strcat(PathName, FileName);
Face = imread(fpath);

Face = imresize(Face, [144 96]);
Face1 = imresize(Face, [288 192]);

axes(handles.axes1) % Select the proper axes
box on
imshow(Face1);

axes(handles.axes2) % Select the proper axes
barh([1 5],'EdgeColor', 'none', 'FaceColor', 'none');
set(handles.axes2, 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], 'Color', [0.9 0.9 0.9]);
box on;

set(handles.bttRecognize, 'Enable', 'on');
set(handles.txtEmotion, 'String', '');

set(handles.txtEmotion, 'String', '');
set(handles.txtProcessTime, 'String', '');


% --- Executes on button press in bttRecognize.
function bttRecognize_Callback(hObject, eventdata, handles)
% hObject    handle to bttRecognize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Face;

individual = load('Best Parameters.mat');
individual = individual.BestIndividual.Parameters;
% individual = [1 5 2 6 1 10 3 15 14 22 9 11 3 11 2 15 2 18 18 24 2 7 3 8 2 15 2 18 18 28 3 5 2 6 7 9];

% file_name = strcat('F:\University Documents\Arshad\Thesis\Face Database\Radboud Emotion\Angry\', int2str(7), '.jpg');
% file_name = strcat('F:\University Documents\Arshad\Thesis\Codes\Training Images\Radboud\New folder\1 (', int2str(i), ').jpg');
% file_name = strcat('F:\University Documents\Arshad\Thesis\Face Database\WSEFEP\Angry\', int2str(i), '.jpg');        
% file_name = strcat('F:\University Documents\Arshad\Thesis\Codes\Test Images\Disgust\', int2str(39), '.jpg');
% file_name = strcat('F:\University Documents\Arshad\Thesis\Codes\Training Images\MMI\', int2str(i), '.jpg');    
% file_name = strcat('F:\University Documents\Arshad\Thesis\Codes\Training Images\Emotions MMI\', int2str(i), '.jpg');    
% file_name = strcat('F:\University Documents\Arshad\Thesis\Codes\Training Images\', int2str(i), '.jpg');    
% file_name = strcat('F:\University Documents\Arshad\Thesis\Codes\Training Images\Correct Features\', int2str(i), '.jpg');    
% file_name = strcat('F:\University Documents\Arshad\Thesis\Codes\Training Images\Training Images\', int2str(70), '.jpg');
% Face = imread(file_name, 'jpg');
% Face = imread('5.jpg', 'jpg');
% Face = imresize(Face, [144 96]);
% figure;
% imshow(Face); 

% color2(Face, 0); color2(Face, 1); color2(Face, 3); color2(Face, 4); color2(Face, 6); color2(Face, 7); color2(Face, 8);  color2(Face, 10);
% Face = color2(Face, 9);

% tic


tic;

%Eye Opening
[EyeOpening1 EyesMidpoint EyesDistance] = Eye_Detection(Face, handles);
% EyeOpening1

% Domain = [0 3; 2 8; 4 12; 10 20; 13 30];
% Parameters = [1 3; 3 8; 2 10; 4 18; 20 22];
Parameters = individual(1: 10);
EyeOpening = Fuzzification(EyeOpening1, Parameters, 5);
% membership_value


%Eyebrow Constriction
[EyebrowConstriction1 EyeWindowUpSide Unimodal] = Eyebrow_Detection(Face, handles);

% if(5*Unimodal+EyebrowConstriction1^2 >= 130)
%     EyebrowSlope = 1;
% else
%     EyebrowSlope = 0;
% end
Parameters = [129 130 129 130];
if(EyebrowConstriction1 ~= -1)
    value = 5*Unimodal+EyebrowConstriction1^2;
else
    value = -1000000;
end
EyebrowSlope = Fuzzification(value, Parameters, 2);

% Domain = [0 10; 8 16; 13 19; 17 23; 18 30];
% Parameters = [8 12; 4 12; 2 18; 1 20; 20 23];
Parameters = individual(11: 20);
% Parameters = [8 12 5 16 5 16 5 16 20 23];
EyebrowConstriction = Fuzzification(EyebrowConstriction1, Parameters, 5);


% Mouth Opening
if(EyeWindowUpSide - 40 > 1)
    Face = Face(EyeWindowUpSide - 40:144, :, :);
else
    Face = Face(1:144, :, :);
end
[MouthOpening1 MouthCornerDisplacement1 UpSide_MouthWindow Mouth_Length Mean_Intensity LipDiameter] = Mouth_Detection(Face, EyesMidpoint, EyesDistance, handles);
% if(Mouth_Length >= 72)
%     MouthLength = 1;
% else
%     MouthLength = 0;
% end

Parameters = [71 72 71 72];
if(MouthOpening1 ~= -1)
    value = Mouth_Length;
else
    value = -1000000;
end
MouthLength = Fuzzification(value, Parameters, 2);

% if(Mean_Intensity >= 325)
%     MeanIntensity = 0;
% else
%     MeanIntensity = 1;
% end
Parameters = [324 325 324 325];
if(MouthOpening1 ~= -1)
    value = Mean_Intensity;
else
    value = -1000000;
end
MeanIntensity = Fuzzification(value, Parameters, 2);

% Domain = [0 7; 4 12; 8 22; 16 26; 18 50];
% Parameters = [2 7; 2 8; 3 15; 2 21; 18 28];
Parameters = individual(21: 30);
MouthOpening = Fuzzification(MouthOpening1, Parameters, 5);
% membership_value

%Mouth Corners
% Domain = [0 7; 4 10; 6 50];
% Parameters = [4 8; 1 8; 7 9];
Parameters = individual(31: 36);
MouthCornerDisplacement = Fuzzification(MouthCornerDisplacement1, Parameters, 3);
% membership_value

%Nose Side Wrinkle
NoseRegion = Face(EyesMidpoint(1) + 10: UpSide_MouthWindow, 20: 75, :);

NoseSideWrinkle = NoseSideWrinkle_Detection(NoseRegion);


% if(5*MouthCornerDisplacement1^2+5*EyeOpening1+3*Mouth_Length+LipDiameter^2-EyebrowConstriction1^2 >= 210)
%     SadBelieve = 1;
% else
%     SadBelieve = 0;
% end
Parameters = [209 210 209 210];
if(MouthOpening1 ~= -1 && EyeOpening1 ~= -1 && EyebrowConstriction1 ~= -1)
    value = 5*MouthCornerDisplacement1^2+5*EyeOpening1+3*Mouth_Length+LipDiameter^2-EyebrowConstriction1^2;
else
    value = -1000000;
end
SadBelieve = Fuzzification(value, Parameters, 2);

[EmotionBelieve Emotion_Intencity Max_Believe] = Fuzzy_Emotion_Recognition(MouthOpening, MouthCornerDisplacement, EyeOpening, EyebrowConstriction, NoseSideWrinkle, MouthLength, MeanIntensity, EyebrowSlope, SadBelieve);
% [EmotionBelieve Emotion_Intencity Max_Believe] = Fuzzy_Emotion_Recognition2(MouthOpening, MouthCornerDisplacement, EyeOpening, EyebrowConstriction, NoseSideWrinkle, MouthLength, MeanIntensity, EyebrowSlope, SadBelieve);

% happy_believe
% sad_believe
% neutral_believe
% fear_believe
% surprise_believe
% angry_believe
% disgust_believeb

% t = toc;

Temp = [MouthOpening1 MouthCornerDisplacement1 EyeOpening1 EyebrowConstriction1 NoseSideWrinkle Mouth_Length Mean_Intensity 5*Unimodal+EyebrowConstriction1^2 5*MouthCornerDisplacement1^2+5*EyeOpening1+3*Mouth_Length+LipDiameter^2-EyebrowConstriction1^2];
% Temp1(i, :) = uint8([happy_believe sad_believe fear_believe surprise_believe angry_believe disgust_believe]*100);

axes(handles.axes2) % Select the proper axes
barh(EmotionBelieve*100, 'FaceColor', [0 0.9 1]);
set(handles.axes2, 'YTickLabel', [], 'TickLength', [0 0], 'Color', [0.9 0.9 0.9]);
box on;

set(handles.txtEmotion, 'String', strcat(Emotion_Intencity, ' :', '  ', num2str(uint8(Max_Believe*100)), '%'));
set(handles.txtProcessTime, 'String', num2str(toc));

clc




%%%%%%%%%%%%%%%%%%%%%%%%%<< Nose Side Wrinkle >>%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NoseSideWrinkle = NoseSideWrinkle_Detection(NoseRegion)

GrayNoseRegion = rgb2gray(NoseRegion);
GrayNoseRegion = medfilt2(GrayNoseRegion);
% figure;
% imshow(GrayNoseRegion);
% GrayNoseRegion = imadjust(GrayNoseRegion);
% figure;
% imshow(GrayNoseRegion);

[row col] = size(GrayNoseRegion);
x = [col/2-10 col/2+11 col/2+11 col/2-10 col/2-10];
y = [1 1 row row 1];
Mask = ~poly2mask(x, y, row, col);
% figure;
% imshow(Mask);

NoseRegionEdges = edge(uint8(GrayNoseRegion), 'log');

NoseRegionEdges = Mask & NoseRegionEdges; 

se = strel([1 1]);
Temp = imopen(NoseRegionEdges, se);
NoseRegionEdges = NoseRegionEdges - Temp;

[NoseRegionEdges num] = bwlabel(NoseRegionEdges, 8);
% figure;
% imshow(NoseRegionEdges);

count = 0;
for i = 1: num
    r = []; c = [];
    [r c] = find(NoseRegionEdges == i);
    [min_r min_ri] = min(r); [max_r max_ri] = max(r); min_c = min(c); max_c = max(c);
    for j = 1: length(max_ri)
        maxr_col(j) = c(max_ri(j));
    end
    for j = 1: length(min_ri)
        minr_col(j) = c(min_ri(j));
    end
    d1 = abs(max(minr_col) - min(maxr_col));
    d2 = abs(min(minr_col) - max(maxr_col));
    width = min(d1, d2);
    
    if((max_r - min_r) < 7 || (max_c - min_c) < 3 || (max_c - min_c) >= 20 || width < 4 || (max_c < col/2 && maxr_col > minr_col) || (min_c > col/2 && maxr_col < minr_col))
        for j = 1: length(r)
            NoseRegionEdges(r(j), c(j)) = 0;
        end
        count = count + 1;
    end
end

% figure;
% subplot(1,2,1);
% imshow(imresize(NoseRegion, 2));
% subplot(1,2,2);
% imshow(NoseRegionEdges);

if(num - count > 0)
    NoseSideWrinkle = 1;
else
    NoseSideWrinkle = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%<< Eye Detection >>%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EyeOpening EyesMidpoint EyesDistance] = Eye_Detection(Face, handles)
global RightEyeCorner;
global LeftEyeCorner;
global RightEyeCenter;
global LeftEyeCenter;

[Row_num Col_num t] = size(Face);
% I = Face(1: Row_num/2, :, :);
% I = rgb2gray(I);
% t = graythresh(I);
% I = im2bw(I, 0.55);
RightEyeRegion1 = Face(Row_num/4: Row_num/2, 1:Col_num/2, :);
LeftEyeRegion1 = Face(Row_num/4: Row_num/2, Col_num/2:Col_num, :);

[RightEyeRegion RightEyeCorner1] = Eye_Region_Detector(RightEyeRegion1);
[LeftEyeRegion LeftEyeCorner1] = Eye_Region_Detector(LeftEyeRegion1);

if(RightEyeRegion == -1)
    RightEyeRegion = RightEyeRegion1;
end
if(LeftEyeRegion == -1)
    LeftEyeRegion = LeftEyeRegion1;
end
% figure;
% subplot(1,2,1);
% imshow(RightEyeRegion);
% subplot(1,2,2);
% imshow(LeftEyeRegion);

[RightEyeCorner LeftEyeCorner] = Eye_Region_Detector_Template_Matching(RightEyeRegion, LeftEyeRegion, RightEyeCorner1, LeftEyeCorner1);

l = 30; w = 15;
end_row = RightEyeCorner(1)+w-1;
if(end_row > size(RightEyeRegion1, 1))
    end_row = size(RightEyeRegion1, 1);
end

end_col = RightEyeCorner(2)+l-1;
if(end_col > size(RightEyeRegion1, 2))
    end_col = size(RightEyeRegion1, 2);
end
RightEyeRegion2 = RightEyeRegion1(RightEyeCorner(1): end_row, RightEyeCorner(2): end_col, :);

end_row = LeftEyeCorner(1)+w-1;
if(end_row > size(LeftEyeRegion1, 1))
    end_row = size(LeftEyeRegion1, 1);     
end

end_col = LeftEyeCorner(2)+l-1;
if(end_col > size(LeftEyeRegion1, 2))
    end_col = size(LeftEyeRegion1, 2);
end
LeftEyeRegion2 = LeftEyeRegion1(LeftEyeCorner(1): end_row, LeftEyeCorner(2): end_col, :);


% I = imcontour(I, 1);
% imshow(I);

% imshow(RightEyeRegion);

% figure;
% subplot(1,2,1);
% imshow(RightEyeRegion2);
% subplot(1,2,2);
% imshow(LeftEyeRegion2);


[row_dis col_dis RightEyeMap] = Eye_Region_Correction(RightEyeRegion2);
RightEyeCorner = ceil(RightEyeCorner + [row_dis(1) col_dis(1)]);

if(RightEyeCorner(1) < 1)
    RightEyeCorner(1) = 1;
end
if(RightEyeCorner(2) < 1)
    RightEyeCorner(2) = 1;
end

[row_dis2 col_dis2 LeftEyeMap] = Eye_Region_Correction(LeftEyeRegion2);
LeftEyeCorner = ceil(LeftEyeCorner + [row_dis2(1) col_dis2(1)]);

if(LeftEyeCorner(1) < 1)
    LeftEyeCorner(1) = 1;
end
if(LeftEyeCorner(2) < 1)
    LeftEyeCorner(2) = 1;
end

end_row = RightEyeCorner(1)+w-1;
if(end_row > size(RightEyeRegion1, 1))
    end_row = size(RightEyeRegion1, 1);
end

end_col = RightEyeCorner(2)+l-1;
if(end_col > size(RightEyeRegion1, 2))
    end_col = size(RightEyeRegion1, 2);
end
RightEyeRegion = RightEyeRegion1(RightEyeCorner(1): end_row, RightEyeCorner(2): end_col, :);

if(get(handles.ShowFacialFeatures, 'Value') == 1)
    axes(handles.axes1); % Select the proper axes
    hold on;
    rectangle('Position',[(RightEyeCorner(2))*2, floor(RightEyeCorner(1)+ Row_num/4-1)*2, (end_col - RightEyeCorner(2))*2, (end_row - RightEyeCorner(1))*2],'Curvature',[0 0],'EdgeColor','y','LineWidth', 2, 'LineStyle', '--');
end

end_row = LeftEyeCorner(1)+w-1;
if(end_row > size(LeftEyeRegion1, 1))
    end_row = size(LeftEyeRegion1, 1);     
end

end_col = LeftEyeCorner(2)+l-1;
if(end_col > size(LeftEyeRegion1, 2))
    end_col = size(LeftEyeRegion1, 2);
end
LeftEyeRegion = LeftEyeRegion1(LeftEyeCorner(1): end_row, LeftEyeCorner(2): end_col, :);

if(get(handles.ShowFacialFeatures, 'Value') == 1)
    axes(handles.axes1); % Select the proper axes
    hold on;
    rectangle('Position',[floor(LeftEyeCorner(2)+Col_num/2)*2, floor(LeftEyeCorner(1)+ Row_num/4 - 1)*2, (end_col - LeftEyeCorner(2))*2, (end_row - LeftEyeCorner(1))*2],'Curvature',[0 0],'EdgeColor','y','LineWidth', 2, 'LineStyle', '--');
end
% figure;
% subplot(1,2,1);
% imshow(RightEyeRegion);
% subplot(1,2,2);
% imshow(LeftEyeRegion);


RightEyeMidpoint_col = RightEyeCorner(2) + (l / 2);
LeftEyeMidpoint_col = LeftEyeCorner(2) + Col_num/2 + (l / 2);
EyesDistance = (LeftEyeMidpoint_col - RightEyeMidpoint_col);
EyesMidpoint = [(Row_num/4 + (RightEyeCorner(1)+ w/2)) (RightEyeMidpoint_col + (EyesDistance/2))];

% figure;
% imshow(Face);
% hold on;
% rectangle('Position',[EyesMidpoint(2) - 2,EyesMidpoint(1) - 2, 2*2, 2*2],'Curvature',[0 0],'EdgeColor','b','LineWidth',2);


% RightEyeRegion = imresize(RightEyeRegion, [30 60]);
% LeftEyeRegion = imresize(LeftEyeRegion, [30 60]);
% figure;
% subplot(1,2,1);
% imshow(RightEyeRegion);
% subplot(1,2,2);
% imshow(LeftEyeRegion);

% Eye_detector(I1);
% Eye_Feature_Points(I1);
% Eye_Feature_Points(LeftEyeRegion);

RightEyeRegion = imresize(RightEyeRegion, [30 60]);
LeftEyeRegion = imresize(LeftEyeRegion, [30 60]);

LeftEyeMap = rgb2ycbcr(LeftEyeRegion);
RightEyeMap = rgb2ycbcr(RightEyeRegion);

[LeftEyeEdges x_left_iris y_left_iris] = Eye_Edges_Detector(LeftEyeRegion, LeftEyeMap);
[RightEyeEdges x_right_iris y_right_iris] = Eye_Edges_Detector(RightEyeRegion, RightEyeMap);

[left_eye_opening left_eye_center]= Eye_Opening(LeftEyeEdges, x_left_iris(1), y_left_iris(1));
[right_eye_opening right_eye_center]= Eye_Opening(RightEyeEdges, x_right_iris(1), y_right_iris(1));

LeftEyeCenter = LeftEyeCorner(1) + left_eye_center;
RightEyeCenter = RightEyeCorner(1) + right_eye_center;

if(left_eye_opening ~= -1)
    if(right_eye_opening ~= -1)
        EyeOpening = (left_eye_opening + right_eye_opening) / 2;
    else
        EyeOpening = left_eye_opening;
    end
else
    if(right_eye_opening ~= -1)
        EyeOpening = right_eye_opening;
    else
        EyeOpening = -1;
    end
end


function [CorrectEyeRegion EyeCorner] = Eye_Region_Detector(EyeRegion)
YCBCR = rgb2ycbcr(EyeRegion);
% [row_num col_num] = size(YCBCR);
% I = EyeRegion(row_num/2: row_num, :, :);
% imshow(I);

% C = makecform('srgb2lab');
% LAB = applycform(EyeRegion, C);
% 
% EyeRegion = rgb2gray(EyeRegion);
% 
% figure;
% imshow(EyeRegion);
% figure;
% imshow(BW_EyeRegion);
se = strel('disk', 2);
Dilate = imdilate(EyeRegion, se);
Erosion = imerode(EyeRegion, se);

[row_num col_num t] = size(EyeRegion);
for i = 1: row_num
    for j = 1: col_num
        EyeMapLum(i, j) = double(Dilate(i, j)) / (double((Erosion(i, j))) + 1);  

        t1 = double(YCBCR(i, j, 2)) .^ 2;
        t2 = double((240 - YCBCR(i, j, 3))) .^ 2;
        t3 = (double(YCBCR(i, j, 2)) ./ double(YCBCR(i, j, 3)));
        EyeMapChr(i, j) = 0.33 .* (t1 + t2 + t3);
        
        EyeMap(i, j) = double(EyeMapLum(i, j)) * double(EyeMapChr(i, j));            
    end
end

Max = max(max(EyeMap));
EyeMap = floor((EyeMap / Max) * 255);
% EyeMap = uint8(EyeMap);
% EyeMap = imadjust(EyeMap);


% figure;
% subplot(1,2,1);
% imshow(EyeRegion);
% subplot(1,2,2);
% imshow(EyeMap, [0 255]);

n = row_num * col_num;
Mean = sum(sum(EyeMap)) / n;
Mean2 = sum(sum(double(EyeMap) .^ 2)) / n;
threshold = (Mean + sqrt(Mean2 - (Mean .^ 2))) / 255;
BwEyeMap = im2bw(uint8(EyeMap), threshold);


se = strel([1 1 1]);
BwEyeMap = imopen(BwEyeMap, se);
se = strel('disk', 1);
BwEyeMap = imdilate(BwEyeMap, se);
% figure;
% imshow(BwEyeMap);

% BwEyeMap = Remove_Small_Area(BwEyeMap, 4);
[BwEyeMap num] = bwlabel(BwEyeMap, 8);

% max_len = 0; label_num = 0;
%     for i = 1: num
%         r = []; c = [];
%         [r c] = find(BwEyeMap == i);
%         len = length(r);
%         if(len > max_len && isempty(find(r == 1, 1)) && isempty(find(r == row_num, 1)) && isempty(find(c == 1, 1)) && isempty(find(c == col_num, 1)))
% %         if(len > max_len)
%             max_len = len;
%             label_num = i;
%         end
%     end

max_val = -100;
% [mr mc] = find(EyeMap == max_val, 1);
if(num > 1)
    label_num = 0;
    for l = 1: num
        r = []; c = [];
        [r c] = find(BwEyeMap == l);        
%         max_val_t = max(max(EyeMap(r, c)));
        max_val_t = sum(sum(EyeMap(r, c)));
%         if(mr <= max(r) && mr >= min(r) && mc <= max(c) && mc >= min(c))
        if( max_val_t > max_val && isempty(find(r <= 2, 1)) && isempty(find(c >= col_num-1, 1)) && isempty(find(c >= 2, 1)))
            label_num = l;
            max_val = max_val_t;
        end
%         Sum = 0;
%         for i = 1: length(r)
%             Sum = Sum + EyeMap(r(i), c(i)); 
%         end
%         values_sum(l) = Sum;
    end
    r = []; c = [];
    [r c] = find(BwEyeMap == label_num);    
    if(label_num == 0 || (~isempty(find(c >= col_num-1, 1)) && ~isempty(find(r <= 2, 1))) || (~isempty(find(c <= 2, 1)) && ~isempty(find(r <= 2, 1))) )
        CorrectEyeRegion = -1;
        EyeCorner = [0 0];
        return;
    end
%     [m, label_num] = max(values_sum);
else
    if(num == 1)
        r = []; c = [];
        [r c] = find(BwEyeMap == 1);    
        if( (~isempty(find(c >= col_num-1, 1)) && ~isempty(find(r <= 2, 1))) || (~isempty(find(c <= 2, 1)) && ~isempty(find(r <= 2, 1))) )
            CorrectEyeRegion = -1;
            EyeCorner = [0 0];
            return;
        end
        
        label_num = 1;
    else
        CorrectEyeRegion = -1;
        EyeCorner = [0 0];
        return;
    end
end

if(label_num > 0)    
    r = []; c = [];
    [r c] = find(BwEyeMap == label_num);
    max_r = max(r);
    min_r = min(r);
    if(min_r - 5 > 0 && max_r + 5 <= row_num)        
%         Temp = EyeRegion(min_r - 5: max_r + 5, 5: col_num, :);
%         [m n] = size(Temp);
        m = 15 - (max_r - min_r + 10);
%         n = 30 - (col_num - 5);
        if(m > 0)
            CorrectEyeRegion = EyeRegion(min_r - 5 - m: max_r + 5, 5: col_num, :);
            EyeCorner = [min_r - 5 - m 5];
        else
            CorrectEyeRegion = EyeRegion(min_r - 5: max_r + 5, 5: col_num, :);
            EyeCorner = [min_r - 5  5];
        end
    else
        CorrectEyeRegion = EyeRegion(min_r: max_r, 5: col_num, :);
        EyeCorner = [min_r 5];
    end
        
else
    CorrectEyeRegion = -1;
    EyeCorner = [0 0];
end

% figure;
% subplot(1,2,1);
% imshow(BwEyeMap);
% subplot(1,2,2);
% imshow(EyeMap, [0 255]);



function [RightEyeCorner LeftEyeCorner] = Eye_Region_Detector_Template_Matching(I1, I2, RightEyeCorner, LeftEyeCorner)

% [Row_num Col_num t] = size(Image);
% I = Image(1: Row_num/2, :, :);
% I = rgb2gray(I);
% % t = graythresh(I);
% % I = im2bw(I, 0.55);
% I1 = I(Row_num/4: Row_num/2, 1:Col_num/2);
% I2 = I(Row_num/4: Row_num/2, Col_num/2:Col_num);

I1 = rgb2gray(I1);
I2 = rgb2gray(I2);

l = 30; w = 15;
% imshow(I1);
template = imread('Eye.jpg', 'jpg');
template = imresize(template, [w l]);
template = rgb2gray(template);

v = -100; vi = 1; vj = 1;
[row_num col_num] = size(I1);
for i = 1: row_num - w
    for j = 1: col_num - l
        Eye_Reg = double(I1(i: i+w-1, j: j+l-1));
%         t = var(var(Eye_Reg));
%         t = std(std(Eye_Reg));
        t = corr2(Eye_Reg, template);
        if(t > v)
            vi = i;
            vj = j;
            v = t;
        end
    end
end

RightEyeCorner = [RightEyeCorner(1)+vi RightEyeCorner(2)+vj];
% Eye_Reg1 = I1(vi: vi+w-1, vj: vj+l-1);

[row_num col_num] = size(I2);
v = -100; vi = 1; vj = 1;
for i = 1: row_num - w
    for j = 1: col_num - l
        Eye_Reg = double(I2(i: i+w-1, j: j+l-1));
%          t = var(var(Eye_Reg));
%         t = std(std(Eye_Reg));
        t = corr2(Eye_Reg, template);
        if(t > v)
            vi = i;
            vj = j;  
            v = t;
        end
    end
end

% Eye_Reg2 = I2(vi: vi+w-1, vj: vj+l-1);

LeftEyeCorner = [LeftEyeCorner(1)+vi LeftEyeCorner(2)+vj];


function [row_dis col_dis YCBCR] = Eye_Region_Correction(EyeRegion)

YCBCR = rgb2ycbcr(EyeRegion);
% [row_num col_num] = size(YCBCR);
% I = EyeRegion(row_num/2: row_num, :, :);
% imshow(I);

% C = makecform('srgb2lab');
% LAB = applycform(EyeRegion, C);
% 
% EyeRegion = rgb2gray(EyeRegion);
% 
% figure;
% imshow(EyeRegion);
% figure;
% imshow(BW_EyeRegion);
se = strel('disk',1);
Dilate = imdilate(EyeRegion, se);
Erosion = imerode(EyeRegion, se);

[row_num col_num t] = size(EyeRegion);
t4 = 0.0; t5 = 0.0;
for i = 1: row_num
    for j = 1: col_num
        EyeMapLum(i, j) = double(Dilate(i, j)) / (double((Erosion(i, j))) + 1);  

        t1 = double(YCBCR(i, j, 2)).^2.0;
        t2 = double((255 - YCBCR(i, j, 3))).^2;
        t3 = (double(YCBCR(i, j, 2))./double(YCBCR(i, j, 3)));
        EyeMapChr(i, j) = 0.33 .* (t1 + t2 + t3);
        
        EyeMap(i, j) = double(EyeMapLum(i, j)) * double(EyeMapChr(i, j));
        
        t4 = t4 + double(YCBCR(i, j, 3)).^2.0;
        t5 = t5 + (double(YCBCR(i, j, 3))./double(YCBCR(i, j, 2))); 
    end
end

% n = 0.95 * (((1/(row_num * col_num)) * t4 / ((1/(row_num * col_num)) * t5)));
% for i = 1: row_num
%     for j = 1: col_num
%         t4 = double(YCBCR(i, j, 3)).^2.0;
%         t5 = (double(YCBCR(i, j, 3))./double(YCBCR(i, j, 2)));
%         MouthMap(i, j) = t4 * (t4 - (n * t5)); 
%     end
% end

se = strel('disk',5);
EyeMap = imdilate(EyeMap, se);

Max = max(max(EyeMap));
EyeMap = floor((EyeMap / Max) * 255);
EyeMap = uint8(EyeMap);
% EyeMap = imadjust(EyeMap);
% EyeMap = im2bw(EyeMap, 0.58);

% figure;
% subplot(1,2,1);
% imshow(EyeRegion);
% subplot(1,2,2);
% imshow(EyeMap, [0 255]);
r = []; c = [];
% [r c] = find(EyeMap > 0.95 * max(max(EyeMap)));
[r c] = find(EyeMap == max(max(EyeMap)));
min_r = min(r);
max_r = max(r);
min_c = min(c);
max_c = max(c);
r_max = (min_r + max_r) / 2;
c_max = (min_c + max_c) / 2;
% EyeMapMax = max(max(EyeMap));
% [r_max c_max] = find(EyeMap == EyeMapMax);

row_dis =  r_max - (row_num/2);
col_dis = c_max - (col_num/2);


function [eye_opening eye_center] = Eye_Opening(EyeEdges, x_iris, y_iris)

EyeEdges = bwlabel(EyeEdges, 8);

i = x_iris;
if(EyeEdges(i, y_iris) == 1)
    i = i - 1;
end
while(EyeEdges(i, y_iris) == 0 && i > 1)
    i = i - 1;
end

if(i >= 1 && EyeEdges(i, y_iris) ~= 0)
%     r = []; c = [];
%     [r c] = find(EyeEdges == EyeEdges(i, y_iris));
%     up_eyelid = min(r);
    c = 1;
    while(EyeEdges(i, y_iris) ~= 0 && i > 1 && c < 3)
        i = i - 1;
        c = c + 1;
    end
    up_eyelid = i; 
    
else
    [L_EyeEdges num] = bwlabel(EyeEdges, 8);
    
    if(num == 0)
        eye_opening = -1;
        eye_center = -1;
        return;
    end 
    
    ind = x_iris;
    while(sum(EyeEdges(ind,:)) == 0 && ind > 1)
        ind = ind - 1;   
    end
    
    count = 0; r_temp = [];
    for i = 1: num
        r = []; c = [];
        [r c] = find(L_EyeEdges == i);
        if(~isempty(find(r == ind, 1)))
            count = count + 1;
            r_temp(count) = min(r);            
        end
    end
    
    up_eyelid = min(r_temp);
         
%     i = x_iris;
%     while(sum(EyeEdges(i,:)) == 0 && i > 1)
%         i = i - 1;   
%     end
%     if(i >= 1 && sum(EyeEdges(i,:)) ~= 0)
% %         r = []; c = [];
% %         c = find(EyeEdges(i,:));
% %         [r c] = find(EyeEdges == EyeEdges(i, c(1)));
% %         up_eyelid = min(r);
%         up_eyelid = i;
%     else
%         eye_opening = -1;
%         eye_center = -1;
%         return;
%     end
end


i = x_iris;
row = size(EyeEdges, 1);
while(EyeEdges(i, y_iris) == 0 && i < row)
    i = i + 1;
end

if(EyeEdges(i, y_iris) ~= 0 && i <= row)
%     r = []; c = [];
%     [r c] = find(EyeEdges == EyeEdges(i, y_iris));
%     down_eyelid = max(r);
    c = 1;
    while(EyeEdges(i, y_iris) ~= 0 && i < row && c < 3)
        i = i + 1;
        c = c + 1;
    end
    down_eyelid = i;
else
    [L_EyeEdges num] = bwlabel(EyeEdges, 8);
    
    if(num == 0)
        eye_opening = -1;
        eye_center = -1;
        return;
    end    
        
    ind = x_iris;
    while(sum(EyeEdges(ind,:)) == 0 && ind < row)
        ind = ind + 1;        
    end
    
    count = 0; r_temp = [];
    for i = 1: num
        r = []; c = [];
        [r c] = find(L_EyeEdges == i);
        if(~isempty(find(r == ind, 1)))
            count = count + 1;
            r_temp(count) = max(r);            
        end
    end
    
    down_eyelid = max(r_temp);
       
    
%     i = x_iris;
%     while(sum(EyeEdges(i,:)) == 0 && i < row)
%         i = i + 1;        
%     end
%     if(sum(EyeEdges(i,:)) ~= 0 && i <= row)
% %         r = []; c = [];
% %         c = find(EyeEdges(i,:));
% %         [r c] = find(EyeEdges == EyeEdges(i, c(1)));
% %         down_eyelid = max(r);
%         down_eyelid = i;
%     else
%         eye_opening = -1;
%         eye_center = -1;
%         return;
%     end
end

eye_opening = down_eyelid - up_eyelid;
eye_center = up_eyelid;% + ceil(eye_opening/2);

        
function [EyeEdges x_iris y_iris] = Eye_Edges_Detector(EyeImage, EyeMap)
% E = EyeImage;
% newmap = brighten(0.5);
% E(:, :) = imadjust(EyeImage(:, :),stretchlim(EyeImage(:, :)),[]);
% % EyeImage1 = imadjust(EyeImage, [0 0 0; 0.35 0.35 0.35], []);
% % EyeImage1 = histeq(EyeImage(:, :), 30);
% % E(:, :) = EyeImage1;
% % figure;
% % subplot(1,2,1);
% % imshow(EyeImage);
% % subplot(1,2,2);
% % imshow(EyeImage1);
% 
% EyeImage = E;

% EyeImage = color2(EyeImage, 10);
% EyeImage = EyeImage + 40;

HSV = rgb2hsv(EyeImage);
% HSV(:, :, 1) = (HSV(:, :, 1) - 0.1);
% EyeImage = hsv2rgb(HSV);
% HSV(:, :, 1) = (HSV(:, :, 1) .* 360) - 180;
% HSV(:, :, 1) = (HSV(:, :, 1) .* 360);
% HSV(:, :, 2) = (HSV(:, :, 2) .* 100);
% HSV(:, :, 3) = (HSV(:, :, 3) .* 100);

HSI = rgb2hsi(EyeImage);

rgb = RGB2rgb(EyeImage);
HSL = rgb2hsl(rgb);
% HSL(:, :, 3) = (HSL(:, :, 3) .* 100);
% imshow(EyeImage);
% EyeMap = imresize(EyeMap, 2);


[m n t] = size(EyeMap);
Sclera1 = zeros(m, n);
Sclera2 = Sclera1; Sclera3 = Sclera1; Sclera4 = Sclera1; Sclera5 = Sclera1; Sclera6 = Sclera1;
for i = 1: m
    for j = 1: n
        if(77 <= EyeMap(i, j, 2) && EyeMap(i, j, 2) <= 127 && 133 <= EyeMap(i, j, 3) && EyeMap(i, j, 3) <= 173)
%         if(EyeMap(i, j, 3) <= 1.5862*EyeMap(i, j, 2)+20 && EyeMap(i, j, 3) >= 0.3448*EyeMap(i, j, 2)+76.2069 && EyeMap(i, j, 3) >= -4.5652*EyeMap(i, j, 2)+234.5652 && EyeMap(i, j, 3) <= -1.15*EyeMap(i, j, 2)+301.75 && EyeMap(i, j, 3) <= -2.2857*EyeMap(i, j, 2)+432.85)
%         if(EyeMap(i, j, 1) > 80 && 85 < EyeMap(i, j, 2) && EyeMap(i, j, 2) < 135 && 135 < EyeMap(i, j, 3) && EyeMap(i, j, 3) < 180)
            Sclera1(i, j) = 1;
        end  
        if(EyeImage(i, j, 1) > 95 && EyeImage(i, j, 1) > 40 && EyeImage(i, j, 1) > 20 && (max([EyeImage(i, j, 1) EyeImage(i, j, 2) EyeImage(i, j, 3)]) - min([EyeImage(i, j, 1) EyeImage(i, j, 2) EyeImage(i, j, 3)])) > 15 && abs(EyeImage(i, j, 1) - EyeImage(i, j, 2)) > 15 && EyeImage(i, j, 1) > EyeImage(i, j, 3) && EyeImage(i, j, 1) > EyeImage(i, j, 3))
            Sclera2(i, j) = 1;
        end
        if(EyeImage(i, j, 1) > 220 && EyeImage(i, j, 1) > 210 && EyeImage(i, j, 1) > 170 && abs(EyeImage(i, j, 1) - EyeImage(i, j, 2)) <= 15 && EyeImage(i, j, 1) > EyeImage(i, j, 3) && EyeImage(i, j, 2) > EyeImage(i, j, 3))
            Sclera2(i, j) = 1;
        end

        if(HSV(i, j, 1) < 0.5)
            HSV(i, j, 1) = HSV(i, j, 1) * 360;
        else
            HSV(i, j, 1) = (HSV(i, j, 1) - 1) * 360;
        end
        
        if((HSV(i, j, 1) >= 0 && HSV(i, j, 1) <= 50) && (HSV(i, j, 2) >= 0.2 && HSV(i, j, 2) <= 0.68) && (HSV(i, j, 3) >= 0.35 && HSV(i, j, 3) <= 1))
            Sclera5(i, j) = 1;
        end
        
        HSV(i, j, 3) = HSV(i, j, 3) * 100;
        if(HSV(i, j, 3) > 40 && HSV(i, j, 2) > 0.2 && HSV(i, j, 2) < 0.6 && (HSV(i, j, 1) < 25 || (HSV(i, j, 1) > 335 && HSV(i, j, 1) < 360)) )
            Sclera6(i, j) = 1;
        end

        HSV(i, j, 2) = HSV(i, j, 2) * 100;
        if(HSV(i, j, 3) >= 40 && HSV(i, j, 1) <= ((-0.4 * HSV(i, j, 3)) + 75) && HSV(i, j, 2) >= 10 && HSV(i, j, 2) <= (-HSV(i, j, 1) - (0.1*HSV(i, j, 3)) + 110) )
            if(HSV(i, j, 1) >= 0)
                if(HSV(i, j, 2) <= ((0.08*(100 - HSV(i, j, 3))*HSV(i, j, 1)) + (0.5 * HSV(i, j, 3))) )
                    Sclera3(i, j) = 1; 
                end
            else
                if(HSV(i, j, 2) <= ((0.5 * HSV(i, j, 1)) + 35) )
                    Sclera3(i, j) = 1; 
                end
            end
        end        
        
%         if(HSI(i, j, 3) > 40)
%             if((HSI(i, j, 2) > 13 && HSI(i, j, 2) < 110) && ((HSI(i, j, 1) > 0 && HSI(i, j, 1) < 28) || (HSI(i, j, 1) > 332 && HSI(i, j, 1) < 360)))
%                 Sclera4(i, j) = 1;
%             end
%             if((HSI(i, j, 2) > 13 && HSI(i, j, 2) < 75) && (HSI(i, j, 1) > 309 && HSI(i, j, 1) < 331))
%                 Sclera4(i, j) = 1;
%             end
%         end
        
%         if(rgb(i, j, 1) / rgb(i, j, 2) > 1.185)
%             sum2 = (rgb(i, j, 1) + rgb(i, j, 2) + rgb(i, j, 3))^2;
%             if(((rgb(i, j, 1)*rgb(i, j, 3)) / sum2) > 0.107 && ((rgb(i, j, 1)*rgb(i, j, 2))/sum2) > 0.112)
%                 Sclera4(i, j) = 1;
%             end
%         end
        if(rgb(i, j, 3) / rgb(i, j, 2) < 1.249)
            sum = (rgb(i, j, 1) + rgb(i, j, 2) + rgb(i, j, 3));
            if(sum / (3*rgb(i, j, 1)) > 0.696 && (1/3 - (rgb(i, j, 3)/sum)) > 0.014 && rgb(i, j, 2)/ (3*sum) < 0.108)
                Sclera4(i, j) = 1;
            end
        end


%         if(HSL(i, j, 1) < 0.5)
%             HSL(i, j, 1) = HSL(i, j, 1) * 360;
%         else
%             HSL(i, j, 1) = (HSL(i, j, 1) - 1) * 360;
%         end
%         if(HSL(i, j, 1) > 20 && HSL(i, j, 1) < 28 && HSL(i, j, 3) < 90)
%             Sclera5(i, j) = 1;
%         end
%         if(abs(EyeImage(i, j, 1) - EyeImage(i, j, 2)) < 10 && abs(EyeImage(i, j, 2) - EyeImage(i, j, 3)) < 10 && abs(EyeImage(i, j, 1) - EyeImage(i, j, 3)) < 10 && EyeImage(i, j, 1) > 0)
%             Sclera5(i, j) = 1;
%         end

    end
end

Big_Sclera = Sclera1 .* Sclera2 .* Sclera3 .* Sclera4 .* Sclera5 .* Sclera6;
Small_Sclera = (Sclera1 .* Sclera2) .* Sclera3;

[r c] = find(Big_Sclera);
if(length(r) < (m*n)/4)
    Sclera = Small_Sclera; 
else
    Sclera = Big_Sclera;
end

Sclera = imfill(~Sclera,'holes'); 

Sclera = Remove_Small_Area(Sclera, 4, 0.33, 8);

se = strel('disk', 1);
% D = imdilate(Sclera, se);
E = imerode(Sclera, se);
% Result = D - E;
Result = Sclera - E;
% se = strel([1 1 1 1]);
% Result = imopen(Result, se);

% figure;
% imshow(EyeMap);

% figure;
% subplot(1,2,1);
% imshow(Sclera4);
% subplot(1,2,2);
% imshow(Result);

Img1 = rgb2gray(EyeImage);

% [n m] = size(Img1);
% Img2 = double(Img1) .^ 2;
% Mean_Img2 = sum(sum(Img2)) / (n*m); 
% Mean_Img1 = sum(sum(Img1)) / (n*m);
% threshold = (Mean_Img1 + (sqrt(Mean_Img2 - (Mean_Img1 .^ 2))) - 10) / 255;
% 
% Sclera = im2bw(Img1, threshold);
% % imshow(Sclera);



Img = medfilt2(Img1);
Img = histeq(Img);

% EyeImg = STD_Based_Mask(Img, 6);
% imshow(EyeImg, [0 1]);


EyeImg = Img1;
Img1 = imerode(Img1, strel([1 1 1 1 1]'));
% r = []; c = [];
% [r c] = find(Img1 <= 100);
% for i = 1: length(r)
%     Img1(r(i), c(i)) = 0;
% end
% figure;
% subplot(1,2,1);
% imshow(EyeImg);
% subplot(1,2,2);
% imshow(Img1);

% se = strel('disk', 1);
% EyeMap = imdilate(EyeMap, se);
% imshow(EyeMap, [0 255]);
% 
% EyeImg = double(EyeMap) - double(Img);
% Max = max(max(EyeImg));
% EyeImg = ceil((EyeImg ./ Max) * 255);
% imshow(EyeImg, [0 255]);

% points = susan(flipud(img'), '-c', '-3');
% draw(img,points,'SUSAN');

% img = rgb2gray(RightEyeRegion);
%     pt  = kp_harris(img);
%     draw(img,pt,'Harris');


% EyeEdges = susan(flipud(img'), '-e', '-3');    
EyeEdges = edge(uint8(Img1), 'canny');
%  EyeEdges = bwmorph(EyeEdges, 'close');

% S = edge(img, 'Sobel', [], 'horizontal');
% % EyeEdges = bwmorph(EyeEdges, 'erode');
% [row_num col_num] = size(img);
% for i = 1: row_num
%     for j = 1: col_num
%         if(S(i, j)==1 && EyeEdges(i, j)==0)
%             T(i, j) = 1;
%         else
%             T(i, j) = EyeEdges(i, j);
%         end
%     end
% end
% % EyeEdges = T;

se = strel([1 1]);
EyeEdges = imopen(EyeEdges, se);
% figure;
% imshow(EyeEdges);

EyeEdges = Remove_Small_Area(EyeEdges, 4, 0.20, 8);

[EyeEdges num] = bwlabel(EyeEdges, 8);
for i = 1: num
    r = []; c = [];
    [r c] = find(EyeEdges == i);
    max_r = max(r); min_r = min(r);
    c2 = [];
    c2 = c(find(r == min_r));
    if(min_r - 1 > 0 && min_r < 15)
%         if(EyeImg(min_r - 1, c2(1)) > 150)
        if(mean(EyeImg(min_r - 1, min(c):max(c))) > 95)    
            for j = 1: length(r)
                EyeEdges(r(j), c(j)) = 0;
            end
        end
    end
    
    c2 = [];
    c2 = c(find(r == max_r));
    if(max_r + 1 < 31 && max_r > 15)
%         if(EyeImg(max_r + 1, c2(1)) > 150)
        if(mean(EyeImg(min_r + 1, min(c):max(c))) > 95)
            for j = 1: length(r)
                EyeEdges(r(j), c(j)) = 0;
            end
        end
    end
end

se = strel([1 1 1 1]);
EyeEdges = imdilate(EyeEdges, se);
% se = strel([1 1]');
% EyeEdges = imopen(EyeEdges, se);

% figure;
% subplot(1,2,1);
% imshow(EyeImg);
% subplot(1,2,2);
% imshow(EyeEdges);


EyeEdges = EyeEdges | logical(Result);

Img(1,1) = 255; Img(1,n) = 255; Img(m,1) = 255; Img(m,n) = 255;

% figure;
% subplot(1,2,1);
% imshow(EyeImage);
% subplot(1,2,2);
% imshow(EyeEdges);
% hold on;

cx = m/2;
cy = n/2;
% [x_iris y_iris] = find(Img == min(min(Img)));
min_val = min(min(Img));
[x_iris y_iris] = find(Img == min_val);
if(length(x_iris) > 1)
    for i = 1: length(x_iris)
        dis(i) = sqrt((x_iris(i) - cx)^2 + (y_iris(i) - cy)^2);
    end
    [min_dis index] = min(dis);
    
    if(min_dis > 5)
        min_val = min_val + 20;
        [x_iris y_iris] = find(Img <= min_val);
        if(length(x_iris) > 1)
            for i = 1: length(x_iris)
                dis(i) = sqrt((x_iris(i) - cx)^2 + (y_iris(i) - cy)^2);
            end
        [min_dis index] = min(dis);
        end
    end
    x_iris = x_iris(index);
    y_iris = y_iris(index);
end
if(EyeEdges(x_iris, y_iris) == 1 && EyeImg(x_iris - 1, y_iris) < EyeImg(x_iris + 1, y_iris))
    x_iris = x_iris - 1;
else
    x_iris = x_iris + 1;
end
% rectangle('Position',[y_iris - 1,x_iris - 1, 2*1, 2*1],'Curvature',[0 0],'EdgeColor','r','LineWidth',2);


function Image = Remove_Small_Area(Image, num, fraction, min_len)

[LabeledImage, n] = bwlabel(Image, num);
length = 0;
for i = 1: n
    [r, c] = find(LabeledImage == i);
    length(i) = max(c) - min(c);
end

length_thershold = max(length) * fraction;

for i = 1: n
    [r, c] = find(LabeledImage == i);
    ln = max(c) - min(c);
    rn = size(r);
    if(ln <= length_thershold || ln <= min_len)
        for k = 1: rn
            Image(r(k), c(k)) = 0;
        end
    end        
end  


function HSI = rgb2hsi(RGB)
[n m t] = size(RGB);
RGB = double(RGB);

for i = 1: n
    for j = 1: m
        I1 = (1/3) * (RGB(i, j, 1) + RGB(i, j, 2) + RGB(i, j, 3));
        I2 = 0.5 * (RGB(i, j, 1) - RGB(i, j, 3));
        I3 = 0.25 *  (2 * RGB(i, j, 2) - RGB(i, j, 1) - RGB(i, j, 3));
        
        temp = atan2(I3, I2);
        if(temp >= 0)
            HSI(i, j, 1) = temp * (180/pi);
        else
            HSI(i, j, 1) = 360 + temp * (180/pi);
        end
        HSI(i, j, 2) = sqrt(I2^2 + I3^2);
        HSI(i, j, 3) = I1;
    end
end


function rgb = RGB2rgb(RGB)
[n m t] = size(RGB);
RGB = double(RGB);

for i = 1: n
    for j = 1: m
        sum = RGB(i, j, 1) + RGB(i, j, 2) + RGB(i, j, 3);
        rgb(i, j, 1) = RGB(i, j, 1) / sum; 
        rgb(i, j, 2) = RGB(i, j, 2) / sum;        
        rgb(i, j, 3) = RGB(i, j, 3) / sum;
    end
end


function HSL = rgb2hsl(rgb)
[n m t] = size(rgb);
rgb = double(rgb);

for i = 1: n
    for j = 1: m
        HSL(i, j, 3) =  rgb(i, j, 1) * 0.29900 + rgb(i, j, 2) * 0.58700 + rgb(i, j, 3) * 0.14400;
        u = - rgb(i, j, 1) * 0.14714 - rgb(i, j, 2) * 0.28886 + rgb(i, j, 3) * 0.43600;
        v =   rgb(i, j, 1) * 0.61500 - rgb(i, j, 2) * 0.51499 - rgb(i, j, 3) * 0.10001;
        HSL(i, j, 1) = atan2( v, u );
        HSL(i, j, 2) = sqrt( u*u + v*v );
    end 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%<< Eyebrow Detection >>%%%%%%%%%%%%%%%%%%%%%%%%
function [EyebrowConstriction EyeWindowUpSide Unimodal] = Eyebrow_Detection(Face, handles)
global RightEyeCorner;
global LeftEyeCorner;
% global RightEyeCenter;
% global LeftEyeCenter;

[Row_num Col_num t] = size(Face);
len = 40; width = 22;
if(RightEyeCorner(1) + Row_num/4 - width > 0 || LeftEyeCorner(1)+ Row_num/4 - width >0)
    RightEyebrowRegion = Face(RightEyeCorner(1)+ Row_num/4 - width: RightEyeCorner(1)+ Row_num/4 +2, 1: Col_num/2, :);
%     LeftEyebrowRegion = Face(LeftEyeCorner(1)+ Row_num/4 -width: LeftEyeCorner(1)+ Row_num/4+2, Col_num/2: Col_num, :);
    LeftEyebrowRegion = Face(RightEyeCorner(1)+ Row_num/4 - width: RightEyeCorner(1)+ Row_num/4 +2, Col_num/2: Col_num, :);
%     RightEyebrowRegion = Face(RightEyeCenter+ Row_num/4 - width: RightEyeCenter+ Row_num/4 -2, 1: Col_num/2, :);
%     LeftEyebrowRegion = Face(LeftEyeCenter+ Row_num/4 -width: LeftEyeCenter+ Row_num/4 -2, Col_num/2: Col_num, :);

    EyeWindowUpSide = RightEyeCorner(1)+ Row_num/4 - width;
else    
    RightEyebrowRegion = Face(1: RightEyeCorner(1)+ Row_num/4+2, 1: Col_num/2, :);
    LeftEyebrowRegion = Face(1: LeftEyeCorner(1)+ Row_num/4+2, Col_num/2: Col_num, :);
%     RightEyebrowRegion = Face(1: RightEyeCenter+ Row_num/4 -2, 1: Col_num/2, :);
%     LeftEyebrowRegion = Face(1: LeftEyeCenter+ Row_num/4 -2, Col_num/2: Col_num, :);

    EyeWindowUpSide = 1;
end

% RightEyebrowRegion = Face(RightEyeCorner(1)-width: RightEyeCorner(1)+5, RightEyeCorner(2)+10: RightEyeCorner(2)+len+5);
% LeftEyebrowRegion = Face(LeftEyeCorner(1)-width: LeftEyeCorner(1)+5, LeftEyeCorner(2)-5: LeftEyeCorner(2)+len-10);

% RightEyebrowRegion = histeq(RightEyebrowRegion);
% LeftEyebrowRegion = histeq(LeftEyebrowRegion);
% RightEyebrowRegion = edge(RightEyebrowRegion, 'sobel');
% LeftEyebrowRegion = edge(LeftEyebrowRegion, 'sobel');

[RightEyebrowConstriction UnimodalR X Y] = Eyebrow_Constriction(RightEyebrowRegion, 'Right');
if(X(1) ~= -1 && Y(1) ~= -1)
    if(get(handles.ShowFacialFeatures, 'Value') == 1)
        axes(handles.axes1); % Select the proper axes
        hold on;
        line(Y*2, (X+EyeWindowUpSide)*2, 'LineWidth',2);
    end
end
[LeftEyebrowConstriction UnimodalL X Y] = Eyebrow_Constriction(LeftEyebrowRegion, 'Left');
if(X(1) ~= -1 && Y(1) ~= -1)
    if(get(handles.ShowFacialFeatures, 'Value') == 1)
        axes(handles.axes1); % Select the proper axes
        hold on;
        line((Y+Col_num/2)*2, (X+EyeWindowUpSide)*2, 'LineWidth',2);
    end
end

Unimodal = UnimodalR + UnimodalL;

if(RightEyebrowConstriction ~= -1 && LeftEyebrowConstriction ~= -1)
    if(abs(RightEyebrowConstriction - LeftEyebrowConstriction) >= 7)
        EyebrowConstriction = min([RightEyebrowConstriction LeftEyebrowConstriction]);
    else
        EyebrowConstriction = ceil((RightEyebrowConstriction + LeftEyebrowConstriction)/2);
    end
else
    if(RightEyebrowConstriction ~= -1)
        EyebrowConstriction = RightEyebrowConstriction;
    else
        if(LeftEyebrowConstriction ~= -1)
            EyebrowConstriction = LeftEyebrowConstriction;
        else
            EyebrowConstriction = -1;
        end
    end
end


function [EyebrowConstriction Unimodal X Y] = Eyebrow_Constriction(EyebrowRegion, EyeSide)

GrayEyebrowRegion = rgb2gray(EyebrowRegion);

len = 40; width = 21;
se = strel([1 1 1 1]');

Dilate = imdilate(GrayEyebrowRegion, se);
Erosion = imerode(GrayEyebrowRegion, se);

MaskImage = Dilate - Erosion;

% Finding Appropriate Threshold
[row col] = size(MaskImage);
s = sum(sum(MaskImage));
Mean = s/(row*col);

M2 = double(MaskImage) .^ 2;
s2 = sum(sum(M2));
Mean2 = s2/(row*col);

threshold = Mean + sqrt(abs(Mean2 - (Mean^2))) + 10;
threshold = (ceil(threshold)) / 255;

% threshold1 = graythresh(MaskImage);
BwMaskImage = im2bw(MaskImage, threshold);
% figure;
% imshow(BwMaskImage);

se = strel([1 1 1]);
BwMaskImage = imdilate(BwMaskImage, se);
% figure;
% imshow(BwMaskImage);

%Remove False Regions
% [BwMaskImage num] = bwlabel(BwMaskImage1, 8);
[BwMaskImage num] = bwlabel(BwMaskImage, 4);
sw = 0; sw1 = 1;
for i = 1: num
    sw1 = 0;
    r = []; c = [];
    [r c] = find(BwMaskImage == i);
%     if((max(c) - min(c)) < len/4)
    max_c = max(c); max_r = max(r); min_c = min(c); min_r = min(r);
    if(((~isempty(find(r < 3, 1)))&& (~isempty(find(c < 3, 1)) || ~isempty(find(c >= col-1, 1)))))
        for j =  1: length(r)
            BwMaskImage(r(j), c(j)) = 0;
        end
        sw1 = 1;
    end 
    if(strcmp(EyeSide, 'Right') && sw1 == 0)
       r1 = floor((max_r - min_r)/2);  
       c1 = floor((max_c - min_c)/2);
       Mean = mean(mean(GrayEyebrowRegion( min_r: min_r+r1, 1: min_c+c1 )));
       if(Mean < 40)
           for j =  1: length(r)
               BwMaskImage(r(j), c(j)) = 0;
           end
       end
       sw = 1;
    end
    if(strcmp(EyeSide, 'Left') && sw1 == 0)
       r1 = floor((max_r - min_r)/2);  
       c1 = floor((max_c - min_c)/2);
       Mean = mean(mean(GrayEyebrowRegion( min_r: min_r+r1, min_c+c1:col )));
       if(Mean < 60)
           for j =  1: length(r)
               BwMaskImage(r(j), c(j)) = 0;
           end
       end
       sw = 1;
    end
end
if(sw == 1)
    [BwMaskImage num] = bwlabel(BwMaskImage, 4);
end
% figure;
% imshow(BwMaskImage);

% [BwMaskImage num] = bwlabel(BwMaskImage, 8);
% [BwMaskImage num] = bwlabel(BwMaskImage, 4);
% RGB = label2rgb(BwMaskImage);
% imshow(BwMaskImage);


count = 1; sw = 0;
Eyebrow_Candidate = struct('Middle_Points', 0, 'min_col', 0, 'max_col', 0);
if(num > 1)
    for i = 1: num
        r = []; c = [];
        [r c] = find(BwMaskImage == i);
        
        counter = 1; 
        Object_Width = 0;
        Middle_Points = 0;
        min_c = min(c); max_c = max(c);
        for j = min_c: max_c
            row_index = find(c == j);
            min_row = r(row_index(1));
            max_row = r(row_index(1));
            for k = 2: length(row_index)
                if(r(row_index(k)) > max_row)
                    max_row = r(row_index(k));
                else
                    if(r(row_index(k)) < min_row)
                        min_row = r(row_index(k));
                    end
                end
            end
            
            Object_Width(counter) = max_row - min_row;
%             Middle_Points(counter) = min_row + ceil(Object_Width(counter)/2);
            Middle_Points(counter) = max_row;
            counter = counter + 1;
        end
        
        object_width = max(Object_Width);
        if(object_width > 2*width/5)
%             for j =  1: length(r)
%                 BwMaskImage(r(j), c(j)) = 0;
%             end
%             sw = 1;
        else
            Eyebrow_Candidate.Middle_Points = Middle_Points;
            Eyebrow_Candidate.min_col = min_c;
            Eyebrow_Candidate.max_col = max_c;
            Eyebrow_Candidates(count) = Eyebrow_Candidate;
            count = count + 1;
        end
%         bw = bwselect(BwMaskImage1, r(ceil(length(r)/2)), c(ceil(length(c)/2)), 8);
%         if((max(r) - min(r)) >= 2*width/5)
%             for j =  1: length(r)
%                 BwMaskImage(r(j), c(j)) = 0;
%             end
%         end
    end
end

if(sw == 1)
    [BwMaskImage num] = bwlabel(BwMaskImage, 8);
%     [BwMaskImage num] = bwlabel(BwMaskImage, 4);
end

% figure;
% imshow(BwMaskImage);
% for i = 1: num
%     r = []; c = [];
%     [r c] = find(BwMaskImage == i);
%     if((max(c) - min(c)) < len/4)
%         for j =  1: length(r)
%             BwMaskImage(r(j), c(j)) = 0;
%         end
%     end    
% end
% [BwMaskImage num] = bwlabel(BwMaskImage, 8);
% figure;
% imshow(BwMaskImage);

if(num == 0)
    EyebrowConstriction = -1;  
    Unimodal = -1;
    X = -1;
    Y = -1;    
    return;
else
    for i = 1: num
        r = []; c = [];
        [r c] = find(BwMaskImage == i);
        ind = find(r == size(BwMaskImage, 1));
        for j = 1: length(ind)
            cols(j) = c(ind(j));
        end
        if(~isempty(ind) && length(r) < 160 && (max(c) - min(c) <= len/2) && ~isempty(find(cols > 10, 1)) && ~isempty(find(cols < 30, 1)))
            for j = 1: length(r)
                BwMaskImage(r(j), c(j)) = 0;
            end            
        end
    end
end
[BwMaskImage num] = bwlabel(BwMaskImage, 8);
% figure;
% imshow(BwMaskImage);

if(num == 0)
    EyebrowConstriction = -1;  
    Unimodal = -1;
    X = -1;
    Y = -1; 
    return;
end

if(num == 1)
    r = []; c = [];
    [r c] = find(BwMaskImage == 1);
        
        counter = 1; 
        Object_Width = 0;
        Middle_Points = 0;
        min_c = min(c); max_c = max(c);
        for j = min_c: max_c
            row_index = find(c == j);
            min_row = r(row_index(1));
            max_row = r(row_index(1));
            for k = 2: length(row_index)
                if(r(row_index(k)) > max_row)
                    max_row = r(row_index(k));
                else
                    if(r(row_index(k)) < min_row)
                        min_row = r(row_index(k));
                    end
                end
            end
            
            Object_Width = max_row - min_row;
            Middle_Points(counter) = min_row + ceil(Object_Width/2);
%             Middle_Points(counter) = max_row;
            counter = counter + 1;
        end
        
            Eyebrow_Candidate.Middle_Points = Middle_Points;
            Eyebrow_Candidate.min_col = min_c;
            Eyebrow_Candidate.max_col = max_c;
            Eyebrow_Candidates(1) = Eyebrow_Candidate; 
            
else
    if(num > 1)
        
        BwMaskImage = Remove_Small_Area(BwMaskImage, 8, 0.4, 8);
        [BwMaskImage num] = bwlabel(BwMaskImage, 8);
        
        if(num == 0)
            EyebrowConstriction = -1;    
            Unimodal = -1;
            X = -1;
            Y = -1; 
            return;
        end
        
        [m n] = size(BwMaskImage);
        min_dis = 100; m = m / 2; n = n / 2;
        for i = 1: num
            r = []; c = [];
            [r c] = find(BwMaskImage == i);
%             state = regionprops(BwMaskImage(min(r):max(r), min(c):max(c)), 'centroid');
%             x = min(r) + ((max(r) - min(r)) / 2);  y = (max(c) - min(c)) / 2; 
%             temp_dis = norm([m-x n-y]);
%             temp_dis = abs(m - x);
            temp_dis = min(r);
            if(temp_dis < min_dis)                
                min_dis = temp_dis;
                min_index = i;
            end                        
        end    
        r = []; c = [];
        [r c] = find(BwMaskImage == min_index);
        min_c = min(c); max_c = max(c);


        counter = 1; Middle_Points = 0;
%         for j = min_c: max_c
%             row_index = find(c == j);
%             min_row = r(row_index(1));
%             for k = 2: length(row_index)
%                 if(r(row_index(k)) < min_row)
%                     min_row = r(row_index(k));
%                 end
%             end
%             
%             Middle_Points(counter) = min_row;
%             counter = counter + 1;             

        for j = min_c: max_c
            row_index = find(c == j);
            min_row = r(row_index(1));
            max_row = r(row_index(1));
            for k = 2: length(row_index)
                if(r(row_index(k)) > max_row)
                    max_row = r(row_index(k));
                else
                    if(r(row_index(k)) < min_row)
                        min_row = r(row_index(k));
                    end
                end
            end
            
            Object_Width = max_row - min_row;
            Middle_Points(counter) = min_row + ceil(Object_Width/2);
%             Middle_Points(counter) = max_row;
            counter = counter + 1;
        end            
            Eyebrow_Candidate.Middle_Points = Middle_Points;
            Eyebrow_Candidate.min_col = min_c;
            Eyebrow_Candidate.max_col = max_c;
            Eyebrow_Candidates(1) = Eyebrow_Candidate;
    end
end


% figure;
% subplot(1,2,1);
% imshow(BwMaskImage);
% subplot(1,2,2);
% imshow(EyebrowRegion);
% hold on;
Y = Eyebrow_Candidates(1).min_col: 1: Eyebrow_Candidates(1).max_col;
X = Eyebrow_Candidates(1).Middle_Points;
% line(Y, X, 'LineWidth',3);  

% EyebrowBottom = RightEyeCorner(1)-width + max(X);
% EyebrowConstriction = RightEyeCenter - EyebrowBottom;
% EyebrowConstriction = max(X);

if(strcmp(EyeSide, 'Right'))
    if(min(X) > size(EyebrowRegion, 1)/2)
        if(length(X)-5 > 0)
            EyebrowConstriction = max(X(length(X)-5: length(X)));
%             if(X(length(X)) - X(length(X)-5) < 3)
%                 Unimodal = 0;
%             else
%                 Unimodal = 1;
%             end            
        else
            EyebrowConstriction = X(length(X));
%             Unimodal = 0;
        end
    else
        EyebrowConstriction = X(length(X));
%         Unimodal = 0;
    end
    Unimodal = X(1) - X(length(X)-2);
else
    if(min(X) > size(EyebrowRegion, 1)/2)
        if(length(X) >= 5)
            EyebrowConstriction = max(X(1: 5));
%             if(X(1) - X(5) < 3)
%                 Unimodal = 0;
%             else
%                 Unimodal = 1;
%             end
        else
            EyebrowConstriction = X(1);
%             Unimodal = 0;
        end        
    else
        EyebrowConstriction = X(1);
%         Unimodal = 0;
    end
   Unimodal = X(length(X)-2) - X(1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%<< Mouth Detection >>%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MouthOpening MouthCornerDisplacement UpSide_MouthWindow Mouth_Length MeanIntensity LipDiameter] = Mouth_Detection(Face, EyesMidpoint, EyesDistance, handles)

row_num1 = size(Face, 1);
MouthRegion = Face(row_num1/2: row_num1, :, :);
% MouthRegion = MouthRegion - 40;
% MouthRegion = color2(MouthRegion, 0);
% MouthRegion = imresize(MouthRegion, 2);
% imshow(MouthRegion);

T = medfilt2(MouthRegion(:,:));
% T = medfilt2(T);
MouthRegion(:,:) = T;
% imshow(MouthRegion);
YCBCR = rgb2ycbcr(MouthRegion);
% HSV = rgb2hsv(MouthRegion);

% map = jet(256);
% newmap = rgb2ycbcr(map);
% imshow(YCBCR(:,:), newmap);


% Temp1 = double(YCBCR(:,:,2)) .^ 2; 
% Max = max(max(Temp1));
% YCBCR(:,:,2) = floor((abs(Temp1) ./ Max) * 255);



[row_num col_num h] = size(MouthRegion);

mid_row = row_num / 2; mid_col = col_num / 2;
t4 = 0.0; t5 = 0.0;
for i = 1: row_num
    for j = 1: col_num        
        t4 = t4 + double(YCBCR(i, j, 3)).^2;
        t5 = t5 + (double(YCBCR(i, j, 3))./double(YCBCR(i, j, 2))); 
    end
end
t4 = (1/(row_num * col_num)) * t4;
t5 = (1/(row_num * col_num)) * t5;
n = 0.95 * (t4 / t5);

% MeanCR = sum(sum(YCBCR(:, :, 3))) / (row_num*col_num);
for i = 1: row_num
    for j = 1: col_num
        t4 = double(YCBCR(i, j, 3)).^2;
        t5 = (double(YCBCR(i, j, 3))./double(YCBCR(i, j, 2)));
        t = t4 * (t4 - (n * t5));        
    
%         distance1 = sqrt(((i - mid_row)/2)^2 + ((j - mid_col)/4)^2) + 1;
%         
%         if(distance1 > 6)
%             MouthMap2(i, j) = t * (1 / (distance1.^0.8));
%         else
%             MouthMap2(i, j) = t * (1 / (distance1.^0.5));
%         end
        
        MouthMap1(i, j) = t;
    end
end
Max1 = max(max(MouthMap1));
MouthMap = floor((abs(MouthMap1) ./ Max1) * 255);
% figure;
% imshow(MouthMap1, [0 Max1]);


Edges = edge(uint8(MouthMap), 'sobel');
% figure;
% imshow(Edges);

se = [1 1 1 1 1]';
Temp = imerode(Edges, se);
Edges = Edges & ~Temp;

se = strel([1 1]);
Edges = imdilate(Edges, se);
Edges = Remove_Small_Area(Edges, 8, 0.2, 10);
% figure;
% imshow(Edges);

% NoseTip_row = ceil((EyesMidpoint(1) + (0.6 * EyesDistance)) - (size(Face, 1)/2));
[Edges num] = bwlabel(Edges, 8);
if(num > 1)
    for i = 1: num
        r = []; c = [];
        [r c] = find(Edges == i);       
    %     if(min(r) < NoseTip_row + 5 || max(r) > row_num - 20)
        if(max(r) > row_num - 20 || max(c) < col_num/2 || min(c) > col_num/2)
            for j = 1: length(r)
                Edges(r(j), c(j)) = 0;
            end
        end
    end
    [Edges num] = bwlabel(Edges, 8);
    if(num > 1)
        for i = 1: num
            r = []; c = [];
            [r c] = find(Edges == i);
            min_rows(i) = min(r);
            max_rows(i) = max(r);
        end
        [v min_ind] = min(min_rows); 
        min_rows = sort(min_rows);        
        if(min_rows(2) - min_rows(1) > 15)
            r = []; c = [];
            [r c] = find(Edges == min_ind);
            if(length(r) < 50)
                for j = 1: length(r)
                    Edges(r(j), c(j)) = 0;
                end
            end
        end        
        [v max_ind] = max(max_rows);
        max_rows = sort(max_rows);
        if(max_rows(num) - max_rows(num - 1) > 8)
            r = []; c = [];
            for k = 1: num
                r = []; c = [];
                [r c] = find(Edges == k);
                if(~isempty(find(r == max_rows(num - 1), 1)))
                    max2_ind = k;
                    break;
                end
            end    
            [r1 c1] = find(Edges == max2_ind);
            [r c] = find(Edges == max_ind);
            if(length(r) < 50 && length(r) < length(r1))
                for j = 1: length(r)
                    Edges(r(j), c(j)) = 0;
                end
            end
        end
    end
end

r = []; c = [];
[r c] = find(Edges);
upper_lip = min(r) - 3;
lower_lip = max(r);
LipDiameter = lower_lip - upper_lip - 3;

if(LipDiameter < 1)
    MouthOpening = -1;
    MouthCornerDisplacement = -1000;
    UpSide_MouthWindow = floor(size(MouthRegion, 1)/2);
    return;
end
% figure;
% imshow(Edges);


% [mr mc] = find(MouthMap1 >= 0.8*Max1);
% for i = 1: size(mr, 1)
%     dis(i) = sqrt(((mr(i) - mid_row))^2 + ((mc(i) - mid_col))^2);
% end
% [v ind] = min(dis);
mid_row = lower_lip + (upper_lip - lower_lip)/2;
for i = 1: row_num
    for j = 1: col_num        
%         distance1 = sqrt(((i - mr(ind))/2)^2 + ((j - mc(ind))/4)^2) + 1;
        distance1 = sqrt(((i - mid_row)/2)^2 + ((j - mid_col)/3)^2) + 1;
        
%         MouthMap2(i, j) = MouthMap1(i, j) * (1 / (distance1));
        if(i > upper_lip && i < lower_lip && j > 16 && j < 80)
            MouthMap2(i, j) = MouthMap1(i, j) * (1 / (distance1.^0.2));
        else
            MouthMap2(i, j) = MouthMap1(i, j) * (1 / (distance1.^0.5));
        end
    end
end

% Mouthmap2 = zeros(row_num, col_num);
% MouthMap2(upper_lip: lower_lip, :) = MouthMap1(upper_lip: lower_lip, :);

Max2 = max(max(MouthMap2));
% figure;
% imshow(MouthMap2, [0 Max2]);

MouthMap2 = floor((abs(MouthMap2) ./ Max2) * 255);

s = sum(sum(MouthMap2));  
Mean = s/(row_num * col_num);
M2 = MouthMap2 .^ 2;
s2 = sum(sum(M2));
Mean2 = s2/(row_num * col_num);
threshold = Mean + 1.6*(sqrt(Mean2 - (Mean^2)));
threshold = (ceil(threshold)) / 255;

% threshold = graythresh(uint8(MouthMap2));

MouthMapT2 = im2bw(uint8(MouthMap2), threshold); 

MouthMapT2 = Remove_Small_Area(MouthMapT2, 4, 0.33, 8);

[MouthMapT2, num] = bwlabel(MouthMapT2);

if(num > 1)
    for i = 1: num
        r = []; c = [];
        [r c] = find(MouthMapT2 == i);
        if(max(r) > 0.85 * row_num)
            for j = 1: length(r)
                MouthMapT2(r(j), c(j)) = 0;
            end
        end
    end
end

% figure;
% imshow(MouthMapT2);
up = 0; down = 0; right = 10; left = 10;

r = []; c = [];
[r c] = find(MouthMapT2);

left_side = min(c) - left;
if(left_side < 1)
    left_side = 1;
end
right_side = max(c) + right;
if(right_side > size(MouthRegion, 2))
    right_side = size(MouthRegion, 2);
end
min_r = min(r);
max_r = max(r);
MouthWindow = MouthRegion(min_r - up: max_r + down, left_side: right_side, :);

if(get(handles.ShowFacialFeatures, 'Value') == 1)
    axes(handles.axes1); % Select the proper axes
    hold on;
    rectangle('Position',[(left_side)*2, (min_r+row_num1/2)*2, (right_side - left_side)*2, (max_r - min_r)*2],'Curvature',[0 0],'EdgeColor','y','LineWidth', 2, 'LineStyle', '--');
end

UpSide_MouthWindow = min_r - up + row_num1/2;

% figure;
% imshow(MouthWindow);
MouthWindow2 = YCBCR(min_r - up: max_r + down, left_side: right_side, :);
MouthMapWindow = MouthMap2(min_r - up: max_r + down, left_side: right_side);
% figure;
% imshow(MouthMapWindow, [0 255]);

s = sum(sum(MouthMapWindow));  
Mean = s/(row_num * col_num);
M2 = MouthMapWindow .^ 2;
s2 = sum(sum(M2));
Mean2 = s2/(row_num * col_num);
threshold = Mean + 2*(sqrt(Mean2 - (Mean^2)));
threshold = (ceil(threshold)) / 255;

BWMouthMap = im2bw(uint8(MouthMapWindow), threshold); 

BWMouthMap = Remove_Small_Area(BWMouthMap, 4, 0.33, 8);

% figure;
% imshow(BWMouthMap);


%%%%%%%%%%%%%%%%%%%%%%%%%%% Mouth Corner Detection %%%%%%%%%%%%%%%%%%%%%%%%
se = strel('disk', 6);
% BwMouth = imdilate(BwMouth, se);
BwMouth = imdilate(BWMouthMap, se);
% figure;
% imshow(BwMouth);

r1 = []; c1 = []; tr = []; tc = [];
[r1 c1] = find(BwMouth);
[tr tc h] = size(BWMouthMap);
MC = zeros(tr, tc);
for i = 1: length(r1)    
%     MC(r(i), c(i)) = double(255 - I1(r(i), c(i), 3)) .^ 2 + double(255 - I1(r(i), c(i), 1)) .^ 4;    
%     MC(r1(i), c1(i)) = double(255 - I1(r1(i), c1(i), 1)) .^ 6;   
    MC(r1(i), c1(i)) = double(255 - MouthWindow2(r1(i), c1(i), 1)) .^ 6;    
end

mc = max(max(MC));
MC = floor((abs(MC) / mc) * 255);

MC = imresize(MC, 2);
BwMouth = imresize(BwMouth, 2);

% figure;
% imshow(MC, [0 255]);

[row col] = size(MC);
MCLeft = MC(:, 1:col/2);
MCRight = MC(:, col/2:col);

BwMouthL = BwMouth(:, 1: col/2);
BwMouthR = BwMouth(:, col/2: col);

% figure;
% imshow(MCRight, [0 255]);
% figure;
% imshow(MCLeft, [0 255]);
% 

%Find Right Mouth Corner
[trr tcr] = find(BwMouthR);

max_val = -100;
for i = 1: length(trr)
    if(MCRight(trr(i), tcr(i)) >= max_val)
       max_val = MCRight(trr(i), tcr(i));
    end
end

n = length(trr);
for i = 1: length(trr)
    t1(n) = trr(i);
    t2(n) = tcr(i);
    n = n - 1;
end
trr = t1; tcr = t2;

i = 1; c = 0; ind = 0;
while(i < length(trr))
    indexes = []; temp = [];
    indexes = find(tcr == tcr(i));
    for j = 1: length(indexes)
        temp(j) = MCRight(trr(indexes(j)), tcr(indexes(j)));
    end

    c = c + 1;
    [max_col_val(c) ind] = max(temp);
    max_col_index(c) = {[trr(indexes(ind)) tcr(indexes(ind))]};
    i = max(indexes) + 1;
end

crx = [];
for i = 1: c
    if(max_col_val(i) >= 0.3 * max_val)
       crx = max_col_index{i}(1);
       cry = max_col_index{i}(2);
       break;
    end
end


%Find Left Mouth Corner
[trl tcl] = find(BwMouthL);

max_val = -100;
for i = 1: length(trl)
    if(MCLeft(trl(i), tcl(i)) >= max_val)
       max_val = MCLeft(trl(i), tcl(i));
    end
end

i = 1; c = 0; max_col_val = []; max_col_index = {[]};
while(i < length(trl))
    indexes = []; temp = [];
    indexes = find(tcl == tcl(i));
    for j = 1: length(indexes)
        temp(j) = MCLeft(trl(indexes(j)), tcl(indexes(j)));
    end

    c = c + 1;
    [max_col_val(c) ind] = max(temp);
    max_col_index(c) = {[trl(indexes(ind)) tcl(indexes(ind))]};
    i = max(indexes) + 1;
end

clx = [];
for i = 1: c
    if(max_col_val(i) >= 0.3 * max_val)
%         mvr = MCRight(trr(i), tcr(i));
       clx = max_col_index{i}(1);
       cly = max_col_index{i}(2);
       break;
    end
end

if(~isempty(clx) && ~isempty(crx))
%     [clx cly] = Mouth_Corner_Correction(MCLeft, clx, cly, MCLeft(clx, cly), 'Left');
%     [crx cry] = Mouth_Corner_Correction(MCRight, crx, cry, MCRight(crx, cry), 'Right');
    cry = cry + col/2;
end

%**************************************************************************
coef = 1.6; 
MouthOpening_prev = 0; MouthOpening_next = 0; count = 1;
while(MouthOpening_prev <= MouthOpening_next && coef < 2.5)
    threshold = Mean + coef*(sqrt(Mean2 - (Mean^2)));
    threshold = (ceil(threshold)) / 255;

    BWMouthMap = im2bw(uint8(MouthMapWindow), threshold); 

    BwMouth = Remove_Small_Area(BWMouthMap, 4, 0.33, 8);

    [BwMouth num] = bwlabel(BwMouth, 8);
    if(num > 1)
        col_num = size(BwMouth, 2);
        max_row = -100;
        for i = 1: num
            r = []; c = [];
            [r c] = find(BwMouth == i);      
            if(max(c) < col_num/2 || min(c) > col_num/2)
                 for j = 1: length(r)
                     BwMouth(r(j), c(j)) = 0;
                 end
            end
            RegionLength(i) = length(r);
            max_r = max(r);
            if(max_r > max_row)
                max_row = max_r;
                index = i;
            end
        end
        for i = 1: num
            if(RegionLength(index) < 2*RegionLength(i)/3)
                r = []; c = [];
                [r c] = find(BwMouth == index);      
                for j = 1: length(r)
                    BwMouth(r(j), c(j)) = 0;
                end
                break;
            end
        end
    end
%     figure;
%     imshow(BwMouth);
%     if(coef == 1.7)
        BwMouth1 = BwMouth;
%     end
    % for i = 1: length(r1)    
    % %     MC(r(i), c(i)) = double(255 - I1(r(i), c(i), 3)) .^ 2 + double(255 - I1(r(i), c(i), 1)) .^ 4;    
    % %     MC(r1(i), c1(i)) = double(255 - I1(r1(i), c1(i), 1)) .^ 6;   
    %     MC(r1(i), c1(i)) = double(255 - MouthWindow2(r1(i), c1(i), 1)) .^ 6;    
    % end
    BwMouth = imresize(BwMouth, 2, 'nearest');
%     figure;
%     imshow(BwMouth);

    Sum = sum(BwMouth, 2);
    % plot(Sum);
    [MouthOpening1 UpperLip1 LowerLip1] = Mouth_Opening1(Sum);
    [MouthOpening2 UpperLip2 LowerLip2 EulerNumber] = Mouth_Opening2(BwMouth);
%     MouthOpening2 = MouthOpening1;
%     UpperLip2 = UpperLip1;
%     LowerLip2 = LowerLip1;
    
%     if(MouthOpening2 == 0 && UpperLip2 == -1000)
    if(EulerNumber == 1)
        MouthOpening(count) = 0;
        UpperLip(count) = UpperLip2;
        LowerLip(count) = LowerLip2;
    else
        if(MouthOpening1 > MouthOpening2)
            MouthOpening(count) = MouthOpening1;
            UpperLip(count) = UpperLip1;
            LowerLip(count) = LowerLip1;
        else
            MouthOpening(count) = MouthOpening2;
            if(MouthOpening1 > 0)
                UpperLip(count) = UpperLip1;
                LowerLip(count) = LowerLip1;            
            else
                UpperLip(count) = UpperLip2;
                LowerLip(count) = LowerLip2;    
            end
        end    

        if(MouthOpening2 > 0 && num == 1 && UpperLip1 > 1 && (Sum(UpperLip1) > Sum(LowerLip1)) && 0.5*MouthOpening1 > MouthOpening2)
            MouthOpening(count) = 0;
        end
    end
    
    MouthOpening_prev = MouthOpening_next;
    MouthOpening_next = MouthOpening(count);
    coef = coef + 0.1;
    count = count + 1;
end
%     figure;
%     imshow(BwMouth);
    
[MouthOpening index] = max(MouthOpening);
UpperLip = UpperLip(index);
LowerLip = LowerLip(index);

if(MouthOpening == 0)
%     se = strel('disk', 3);
%     BwMouth = imdilate(BwMouth, se);
%     figure;
%     imshow(BwMouth);
    Sum = sum(BwMouth, 2);
    MouthCenterLine = find(Sum == max(Sum));        
%     MouthCenterLine = size(MouthMapWindow, 1);
    
    % Mouth Center Line Correction
    n = 5;
    if(2*n+1 > size(MC, 1))
        num = size(MC, 1);
    else
        num = 2*n+1;
    end
    for i = 1: num 
        l = (MouthCenterLine(1)-(n+1) + i);
        if(l <= 0)
            l = 1;
        end
        if(l > size(MC, 1))
            break;
        end
        illuminance_sum(i) = sum(MC(l, :)); 
    end
    [v index] = max(illuminance_sum);
    MouthCenterLine = ((MouthCenterLine(1)-(n+1)) + index);  
    if(MouthCenterLine > 3*size(MC, 1)/4)
        MouthCenterLine = size(MC, 1)/2;
    end
    MeanIntensity = -1;
    
else
    if(MouthOpening == -1)
        MouthCornerDisplacement = -1;
        return;
    end
        
    MouthCenterLine = (floor((UpperLip + LowerLip) / 2));
        % Mouth Center Line Correction
    n = 5;
    if(2*n+1 > size(MC, 1))
        num = size(MC, 1);
    else
        num = 2*n+1;
    end
    for i = 1: num 
        l = (MouthCenterLine(1)-(n+1) + i);
        if(l <= 0)
            l = 1;
        end
        if(l > size(MC, 1))
            break;
        end
        illuminance_sum(i) = sum(MC(l, :)); 
    end
    [v index] = max(illuminance_sum);
    MouthCenterLine2 = ((MouthCenterLine(1)-(n+1)) + index);
    
    Gray_Mouth = rgb2gray(MouthWindow);
%     Gray_Mouth = imadjust(Gray_Mouth);
    c1 = floor(cly/2+(cry -  cly)/4)-7; c2 = floor(cly/2+(cry -  cly)/4)+7;
    Mouth_Box = Gray_Mouth(UpperLip/2: LowerLip/2, c1: c2);
    BwMouth_Box = BwMouth1(UpperLip/2: LowerLip/2, c1: c2);
       
    r = []; c = [];
    [r c] = find(BwMouth_Box == 0);
    r2 = []; c2 = [];
    [r2 c2] = find(BwMouth_Box == 1);
%     MeanIntesity = mean(mean(Mouth_Box));
%     MeanIntensity = length(find(Mouth_Box < 50));
    [BwMouth, num] = bwlabel(BwMouth1);
    Temp = regionprops(BwMouth1, 'EulerNumber');        
    if(num == 1 && Temp.EulerNumber == 1 && sum(BwMouth_Box(1,:)) < size(BwMouth_Box, 2)/2)        
        count = length(find(BwMouth_Box));
    else
        count = 0;
        for i = 1: length(r)
            if(Mouth_Box(r(i), c(i)) < 70)
                count = count + 1;
            end
        end
    end
    
    max_mc = max(max(MC));
    mean_mc = mean(mean(MC));
    c1 = c1*2; c2 = c2*2;
    MeanIntensity = length(find(MC(1: size(MC, 1), c1: c2) > 0.5*max_mc));
%     MeanIntensity = length(find(MC(1: size(MC, 1), c1: c2) > mean_mc));
    MeanIntensity = MeanIntensity + count^2;
    if(MeanIntensity >= 52)
        MouthCenterLine = MouthCenterLine2;
    end
    
%     figure;
%     subplot(1,2,1);
%     imshow(MC(1: size(MC, 1), c1: c2), [0 255]);
%     subplot(1,2,2);
%     imshow(BwMouth1);
end

if(~isempty(clx) && ~isempty(crx))
    if(clx < MouthCenterLine(1) && crx < MouthCenterLine(1))
    %     MouthCornerDisplacement = MouthCenterLine(1) - min(clx, crx);
        MouthCornerDisplacement = MouthCenterLine(1) - ((clx + crx)/2);
    else
        if(clx > MouthCenterLine(1) && crx > MouthCenterLine(1))
    %         MouthCornerDisplacement = MouthCenterLine(1) - max(clx, crx);
            MouthCornerDisplacement = MouthCenterLine(1) - ((clx + crx)/2);
        else
            tx1 = MouthCenterLine(1) - clx;
            tx2 = MouthCenterLine(1) - crx;        
            MouthCornerDisplacement = (tx1 + tx2) / 2;
        end
    end
else
    if(isempty(crx) && isempty(crx))
        MouthCornerDisplacement = -1;
    else
        if(isempty(crx))
            MouthCornerDisplacement = clx;
        else
            MouthCornerDisplacement = crx;
        end
    end
end

Mouth_Length = cry -  cly;

MouthWindow = imresize(MouthWindow, 2);

% figure;
% subplot(1,2,1);
% imshow(MouthRegion);
% imshow(MouthWindow);
% hold on;
% s = 1.2;
% rectangle('Position',[cly(1)-s, clx(1)-s, 2*s, 2*s],'Curvature',[0 0],'EdgeColor','y','LineWidth',2);
% rectangle('Position',[cry(1)-s, crx(1)-s, 2*s, 2*s],'Curvature',[0 0],'EdgeColor','y','LineWidth',2);
% 
% Y = [MouthCenterLine(1), MouthCenterLine(1)]; X = [1, size(MouthWindow, 2)];
% line(X, Y, 'Color','y');
% 
% subplot(1,2,2);
% imshow(MC, [0 255]);
% hold on;
% s = 1.2;
% rectangle('Position',[cly(1)-s, clx(1)-s, 2*s, 2*s],'Curvature',[0 0],'EdgeColor','y','LineWidth',2);
% rectangle('Position',[cry(1)-s, crx(1)-s, 2*s, 2*s],'Curvature',[0 0],'EdgeColor','y','LineWidth',2);
% 
% Y = [MouthCenterLine(1), MouthCenterLine(1)]; X = [1, size(MouthWindow, 2)];
% line(X, Y, 'Color','y');



% Mouth_Corrner_Detector(MouthRegion, MouthMap);
 
% figure;
% subplot(1,2,1);
% imshow(MouthRegion);
% subplot(1,2,2);
% imshow(BwMouth, [0 1]);


function [MouthOpening UpperLip LowerLip] = Mouth_Opening1(Sum)

count = 1; count1 = 0; count2 = 0; s = 2;
row = length(Sum);
Sum(row + 1) = 0; 
row = row + 1;
% figure
% plot(Sum);
% LocalMin(count2) = {[1, 0]};

for i = 2: row - 1
    Sum(i) = floor((Sum(i-1) + Sum(i+1)) / 2);
%     Sum(i) = floor((Sum(i-2) + Sum(i-1) + Sum(i+1) + Sum(i+2)) / 4);    
end

% figure;
% plot(Sum);

max_sum = max(Sum); sw = 0;

% Sum = Sum - min(Sum);

for i = 1: row
    if(Sum(i) < 0.25 * max_sum)
        Sum(i) = 0;
    end
end

count2 = count2 + 1;
LocalMin(count2) = {[1, 0]};

% if(Sum(2) > 0.5*max_sum)            
%             count1 = count1 + 1;
%             LocalMax(count1) = {[2, Sum(2)]};
%             s = 3;
% end

if(Sum(2) >= 0)
    i = 2;
    while(~(Sum(i) < Sum(i-1) && Sum(i) <= Sum(i+1)))
        if((Sum(i) >= Sum(i-1) && Sum(i) > Sum(i+1)))
            sw = 1;
            count1 = count1 + 1;
            LocalMax(count1) = {[i Sum(i)]};
            s = i + 1;
            break;
        end
        i = i + 1;
        if(i >= length(Sum))
            break;
        end
    end
    if(sw == 0)
        if(Sum(2) > 0.4*max_sum)
            count1 = count1 + 1;
            LocalMax(count1) = {[2, Sum(2)]};        
        else
            s = i + 1;
        end        
    end
end

% if(Sum(2) >= 0)
%     i = 2;
%     while(Sum(i) ~= 0)
%         if(Sum(i) >= Sum(i-1) && Sum(i) > Sum(i+1))
%             sw = 1;
%             count1 = count1 + 1;
%             LocalMax(count1) = {[i Sum(i)]};
%             s = i + 1;
%             break;
%         end
%         i = i + 1;        
%     end
%     if(sw == 0)
%         if(Sum(2) > 0.5*max_sum)
%             count1 = count1 + 1;
%             LocalMax(count1) = {[2, Sum(2)]};        
%         end
%         i = 2;
%         while(Sum(i) > Sum(i+1))            
%             i = i + 1;
%         end
%         
%         count2 = count2 + 1;
%         LocalMin(count2) = {[i, 0]};
%         s = i + 1;
%     else
%         count2 = count2 + 1;
%         LocalMin(count2) = {[1, 0]};
%         
% %         i = 2;
% %         while(i < row - 2)
% %             if(Sum(i) < Sum(i-1) && Sum(i) <= Sum(i+1))              
% %                 count2 = count2 + 1;
% %                 LocalMin(count2) = {[i, Sum(i)]};
% %                 s = i + 1;
% %                 break;
% %             end
% %             i = i + 1;        
% %         end        
%     end
% else    
%     count2 = count2 + 1;
%     LocalMin(count2) = {[1, 0]};
%     s = 3;
% %     if(Sum(5) == 0)
% %         count2 = count2 + 1;
% %         LocalMin(count2) = {[5, 0]};
% %         s = 6;
% %     else
% %         i = 6;
% %         while(Sum(i) > Sum(i))            
% %             i = i + 1;
% %         end
% %         
% %         count2 = count2 + 1;
% %         LocalMin(count2) = {[i, 0]};
% %         s = i + 1;
% %     end
% end

% [local_min index] = findpeaks(-Sum,'minpeakdistance',2);
% for i = 1: length(local_min)
%     count2 = count2 + 1; 
%     LocalMin(count2) = {[index(i) local_min(i)]};
% end

for i = s: row - 1
    if(Sum(i) >= Sum(i-1) && Sum(i) > Sum(i+1))
        count1 = count1 + 1;
        LocalMax(count1) = {[i Sum(i)]};    
    else
        if(Sum(i) < Sum(i-1) && Sum(i) <= Sum(i+1))
%             if(count1 > 0 && Sum(i) <= (2*LocalMax{count1}(2)/3))                
%                 count2 = count2 + 1;
%                 LocalMin(count2) = {[i Sum(i)]};
%             end
            count2 = count2 + 1;
            LocalMin(count2) = {[i Sum(i)]};
        end
    end
end
  
% l = length(LocalMax);
% for i = 1: l
%     if(LocalMax{i}(2) < 0.25 * max_sum)
%         for j = i: l - 1
%             LocalMax(j) = LocalMax(j+1);
%         end
%         LocalMax(l) = [];
%     end
% end



% Temp(count) = {[1 0]};
% for i = 1: length(LocalMin)-1
%     if(LocalMin{i}(2) <= (2*LocalMax{i}(2)/3) && LocalMin{i}(2) < (2*LocalMax{i+1}(2)/3))        
%         count = count + 1;
%         Temp(count) = LocalMin(i);
%     end
% end
% if(length(LocalMin) == 1)
%     count = count + 1;
%     Temp(count) = LocalMin(1);
% end
% 
% LocalMin = Temp;

% figure;
% plot(Sum);

if(length(LocalMax) == length(LocalMin))
    LocalMin(length(LocalMin)+1) = {[row 0]};
end

LocalMin_T(1) = LocalMin(1);
count = 1;
for i = 1: length(LocalMax)
    if(LocalMax{i}(2) > 0.3*max_sum)
        LocalMax_T(count) = LocalMax(i);
        LocalMin_T(count+1) = LocalMin(i+1);
        count = count + 1;
    end
end
LocalMax = LocalMax_T;
LocalMin = LocalMin_T;

if(length(LocalMax) >= length(LocalMin))
    UpperLip = -1;
    LowerLip = -1;        
    MouthOpening = -1;
    return;
end

False_Max_Min = struct('index', 0, 'orientation', 'right');

c = 0; Max_Temp = False_Max_Min;
for i = 1: length(LocalMax)
    if(LocalMin{i+1}(2) >= 0.75*LocalMax{i}(2))
%     if((LocalMax{i}(2) - LocalMin{i}(2)) >= 0.4*(LocalMax{i}(2) - LocalMin{i+1}(2)) || (LocalMin{i+1}(2) > 0.7*LocalMax{i}(2)))
        c = c + 1;
        False_Max_Min.index = i;
        False_Max_Min.orientation =  'right';
        Max_Temp(c) = False_Max_Min;
        continue;        
    end
    if(LocalMin{i}(2) >= 0.75*LocalMax{i}(2))
%     if((LocalMax{i}(2) - LocalMin{i}(2)) <= 0.4*(LocalMax{i}(2) - LocalMin{i+1}(2)) || (LocalMin{i}(2) > 0.7*LocalMax{i}(2)))
        c = c + 1;
        False_Max_Min.index = i;
        False_Max_Min.orientation =  'left';
        Max_Temp(c) = False_Max_Min;
    end
end

False_Max_Min.index = 0;
c = 0; Min_Temp = False_Max_Min;
for i = 2: length(LocalMin) - 1
    if(LocalMin{i}(2) >= 0.75*LocalMax{i}(2))
%     if(LocalMax{i}(2) < 0.3*LocalMax{i-1}(2))
%     if((LocalMax{i}(2) - LocalMin{i}(2)) <= 0.4*(LocalMax{i-1}(2) - LocalMin{i}(2)) || (LocalMax{i}(2) < 0.3*LocalMax{i-1}(2)))
        c = c + 1;
        False_Max_Min.index = i;
        False_Max_Min.orientation =  'right';
        Min_Temp(c) = False_Max_Min;
        continue;
    end
    if(LocalMin{i}(2) >= 0.75*LocalMax{i-1}(2))
%     if(LocalMax{i-1}(2) < 0.3*LocalMax{i}(2))
%     if((LocalMax{i-1}(2) - LocalMin{i}(2)) <= 0.4*(LocalMax{i}(2) - LocalMin{i}(2)) || (LocalMax{i-1}(2) < 0.3*LocalMax{i}(2)))
        c = c + 1;
        False_Max_Min.index = i;
        False_Max_Min.orientation =  'left';
        Min_Temp(c) = False_Max_Min;
    end
end
False_Max_Min.index = 0;
Min_Temp = False_Max_Min;


lm_sw = 0;
if(Max_Temp(1).index ~= 0)
    c = 0; temp_index = [];
    for i = 1: length(Max_Temp)
        if(strcmp(Max_Temp(i).orientation, 'right'))
            for j = 1: length(LocalMin)
                if(LocalMin{j}(1) > LocalMax{Max_Temp(i).index}(1))
                    c = c + 1;
                    temp_index(c) = j;
                    break;
                end
            end
        else
            for j = length(LocalMin) - 1: -1: 1
                if(LocalMin{j}(1) < LocalMax{Max_Temp(i).index}(1))
                    c = c + 1;
                    temp_index(c) = j;
                    break;
                end
            end            
        end
    end

    if(~isempty(temp_index))
        c = 0; lm_sw = 1;        
        for i = 1: length(LocalMin)
            if(isempty(find(temp_index == i, 1)))
                c = c + 1;
                LocalMinTemp(c) = LocalMin(i);
            end
        end
%         LocalMin = Temp;
    end
end

if(Min_Temp(1).index ~= 0)
    c = 0; temp_index = [];
    for i = 1: length(Min_Temp)
        if(strcmp(Min_Temp(i).orientation, 'right'))
            for j = 1: length(LocalMax)
                if(LocalMax{j}(1) > LocalMin{Min_Temp(i).index}(1))
                    c = c + 1;
                    temp_index(c) = j;
                    break;
                end
            end
        else
            for j = length(LocalMax): -1: 1
                if(LocalMax{j}(1) < LocalMin{Min_Temp(i).index}(1))
                    c = c + 1;
                    temp_index(c) = j;
                    break;
                end
            end
        end
    end
    
    if(~isempty(temp_index))
        c = 0;
        for i = 1: length(LocalMax)
            if(isempty(find(temp_index == i, 1)))
                c = c + 1;
                Temp(c) = LocalMax(i);
            end
        end
        LocalMax = Temp;
    end
end

if(lm_sw == 1)
    LocalMin = LocalMinTemp;
end


LM = {[3 Sum(3)]}; 
for i = 1: length(LocalMin)-1
    Temp = {[0 0]}; count = 0;
    for j = 1: length(LocalMax)
        if(LocalMax{j}(1) > LocalMin{i}(1) && LocalMax{j}(1) < LocalMin{i+1}(1))
            count = count + 1;
            Temp(count) = LocalMax(j); 
        end
    end
    
    if(Temp{1}(1) ~= 0)
        Max = {[1 -1]};
        for k = 1: length(Temp)
            if(Temp{k}(2) >= Max{1}(2))
                Max = Temp(k);
            end
        end
        LM(i) = Max;
    end
end

LocalMax = LM;

c = 0;
for i = 1: length(LocalMin)-1
    if(LocalMin{i}(2) <= 0.75*LocalMax{i}(2))
        c = c + 1;
        Temp2(c) = LocalMax(i);
    end
end
if(c > 0)
    LocalMax = Temp2;
end

if(length(LocalMax) >= 3)
    i = LocalMax{3}(1);
    if(Sum(i) == Sum(i-1))
        while(Sum(i) == Sum(i-1))
            i = i - 1;
        end
    end
    LowerLip = i;
    
    j = LocalMax{1}(1);
%     while(Sum(j) >= 0.98 * LocalMax{1}(2))
%         j = j + 1;
%     end
    UpperLip = j;
    
    MouthOpening = LowerLip - UpperLip;
else
    if(length(LocalMax) == 2)
        i = LocalMax{2}(1);
        if(Sum(i) == Sum(i-1))
            while(Sum(i) == Sum(i-1))
                i = i - 1;
            end
        end
        LowerLip = i;
        
        j = LocalMax{1}(1);
%         while(Sum(j) >= 0.98 * LocalMax{1}(2))
%             j = j + 1;
%         end
        UpperLip = j;
        MouthOpening = LowerLip - UpperLip;
    else
        UpperLip = 0;
        LowerLip = 0;        
        MouthOpening = 0;
    end    
end


function [MouthOpening UpperLip LowerLip EulerNumber] = Mouth_Opening2(BwMouth)
BwMouth = ~Remove_Small_Area(~BwMouth, 4, 0.0001, 5);

BwMouth = padarray(BwMouth, [1 1]);

% num_line = 40;
[BwMouth num] = bwlabel(BwMouth, 8);
euler_num = regionprops(BwMouth, 'EulerNumber');
Boundaries = bwboundaries(BwMouth, 8, 'noholes');

EulerNumber = 0;
if(num == 1 && euler_num.EulerNumber == 1 && length(Boundaries{1}) < 400)
%     MouthOpening = 0;
%     UpperLip = -1000;
%     LowerLip = -1000;
%     return;
    EulerNumber = 1;
end
sw1 = 0;
if(num == 1 && euler_num.EulerNumber < 1)
    sw1 = 1;
end

if(num == 1)
    [r c] = find(BwMouth);
    min_col = min(c);
    min_row = min(r);
    max_c = max(c);
    num_line = (max_c - min_col);
%     interval = floor((max_c - min_col)/num_line);
    interval = 1; 
else
    min_row = 1000;
    for i = 1: num
        r = []; c = [];
        [r c] = find(BwMouth == i);
        min_r = min(r);
        if(min_r < min_row)
            min_row = min_r;
            component_num = BwMouth(r(1), c(1));            
        end
    end        
    r = []; c = [];
    [r c] = find(BwMouth == component_num);        
    min_col = min(c);
    max_c = max(c);
    num_line = (max_c - min_col);
%     interval = floor((max(c) - min_col)/num_line);
    interval = 1;
end

for i = 1: num_line    
    col = min_col + (i-1)*(interval);
    row = min_row;
    while(BwMouth(row, col) == 0)
        row = row + 1;
    end
    while(BwMouth(row, col) > 0)
        row = row + 1;
    end
    Upper_Lip(i) = row - 1;
    row = size(BwMouth, 1) - 1;
    while(row > 0 && BwMouth(row, col) == 0)
        row = row - 1;
    end
    while(row > 0 && BwMouth(row, col) > 0)
        row = row - 1;        
    end
    Lower_Lip(i) = row + 1;
    
    if(num == 1 && length(Boundaries{1}) < 400)
        sw = 0;
        for j = 1: length(Boundaries{1})
            if((Boundaries{1}(j, 1) == Upper_Lip(i) && Boundaries{1}(j, 2) == col) || (Boundaries{1}(j, 1) == Lower_Lip(i) && Boundaries{1}(j, 2) == col))
                sw = 1;
                break;
            end
        end
        if(sw == 1)
            Mouth_Opening(i) = -1000;
        else
            if(Lower_Lip(i) < size(BwMouth, 1)/2)
                Mouth_Opening(i) = 0;
                Upper_Lip(i) = -1;
                Lower_Lip(i) = -1;                
            else
                Mouth_Opening(i) = Lower_Lip(i) - Upper_Lip(i);     
            end
        end        
    else
        if((Upper_Lip(i) < size(BwMouth, 1)/3 || Lower_Lip(i) > 2*size(BwMouth, 1)/3) && length(Boundaries{1}) < 400)
            sw = 0;
            for j = 1: length(Boundaries{1})
                if((Boundaries{1}(j, 1) == Upper_Lip(i) && Boundaries{1}(j, 2) == col) || (Boundaries{1}(j, 1) == Lower_Lip(i) && Boundaries{1}(j, 2) == col))
                    sw = 1;
                    break;
                end
            end
            if(sw == 1)
                Mouth_Opening(i) = -1000;
            else            
                Mouth_Opening(i) = Lower_Lip(i) - Upper_Lip(i);
            end
        else
            Mouth_Opening(i) = Lower_Lip(i) - Upper_Lip(i);
        end
    end
end

num1 = find(Lower_Lip == -1);
num2 = find(Mouth_Opening == -1000);
if(sw1 == 1 && length(num1) + length(num2) == num_line)
    [BwMouth2 num] = bwlabel(~BwMouth, 8); 
    for i = 1: num        
        r = []; c = [];
        [r c] = find(BwMouth2 == i);
        if(max(r) < size(BwMouth2, 1)/2 && length(r) < 70 && (max(c) - min(c) < 15))
            EulerNumber = 1;
        end
    end    
end

[MouthOpening index] = max(Mouth_Opening);
UpperLip = Upper_Lip(index);
LowerLip = Lower_Lip(index);
col = min_col + (index-1)*(interval);

% figure;
% imshow(BwMouth);
% hold on;
% rectangle('Position',[col - 1,UpperLip - 1, 2*1, 2*1],'Curvature',[0 0],'EdgeColor','b','LineWidth',2);
% rectangle('Position',[col - 1,LowerLip - 1, 2*1, 2*1],'Curvature',[0 0],'EdgeColor','r','LineWidth',2);



    
    
    
%%%%%%%%%%%%%%%%%%%%%<< Fuzzy Emotion Recognition >>%%%%%%%%%%%%%%%%%%%%%%%%
function membership_value = Fuzzification(value, Parameters, mf_num)

switch(mf_num)
    case 2
        if(value == -1000000)
            membership_value = [1 1];
%             sign = 'unsigned';
            return;
        end
        
        membership_value(1) = evalmf(value, [Parameters(1) Parameters(2)], 'zmf');
        membership_value(2) = evalmf(value, [Parameters(3) Parameters(4)], 'smf');
        return;
        
    case 3
        if(value == -1000)
            membership_value = [1 1 1 1 1 1];
%             sign = 'unsigned';
            return;
        end
        
        if(value < 0)
%             sign = 'Negative';
            value = value * -1;            
            index = 3;
        else
            index = 0;
%             sign = 'Positive';
        end
%         value = floor(value);
        membership_value = [0 0 0 0 0 0];
        
%         x = Domain(1, 1):1:Domain(1, 2);
%         Very_Small = zmf(x, [Parameters(1, 1) Parameters(1, 2)]);
% %         figure;
% %         plot(x, Very_Small);
%         if(value >= Domain(1, 1) && value <= Domain(1, 2))
%             membership_value(index + 1) = Very_Small(value - (Domain(1, 1)-1));
%         else
%             membership_value(index + 1) = 0;
%         end
% 
%         x = Domain(2, 1):1:Domain(2, 2);
%         Small = gaussmf(x,[Parameters(2, 1) Parameters(2, 2)]);
% %         hold on;
% %         plot(x, Small);
%         if(value >= Domain(2, 1) && value <= Domain(2, 2))
%             membership_value(index + 2) = Small(value - (Domain(2, 1)-1));
%         else
%             membership_value(index + 2) = 0;
%         end
% 
%         x = Domain(3, 1):1:Domain(3, 2);
%         Medium = smf(x,[Parameters(3, 1) Parameters(3, 2)]);
% %         hold on;
% %         plot(x, Medium);
%         if(value >= Domain(3, 1) && value <= Domain(3, 2))
%             membership_value(index + 3) = Medium(value - (Domain(3, 1)-1));
%         else
%             membership_value(index + 3) = 0;
%         end

        membership_value(index + 1) = evalmf(value, [Parameters(1) Parameters(2)], 'zmf');
        membership_value(index + 2) = evalmf(value, [Parameters(3) Parameters(4)], 'gaussmf');
        membership_value(index + 3) = evalmf(value, [Parameters(5) Parameters(6)], 'smf');
        return;
        
    case 5   
        if(value == -1)
            membership_value = [1 1 1 1 1];            
            return;
        end
              
%         value = floor(value);
        
%         x = Domain(1, 1):1:Domain(1, 2);
%         Very_Small = zmf(x, [Parameters(1, 1) Parameters(1, 2)]);
% %         figure;
% %         plot(x, Very_Small);
%         if(value >= Domain(1, 1) && value <= Domain(1, 2))
%             membership_value(1, 1) = Very_Small(value - (Domain(1, 1)-1));
%         else
%             membership_value(1, 1) = 0;
%         end
% 
%         x = Domain(2, 1):1:Domain(2, 2);
%         Small = gaussmf(x,[Parameters(2, 1) Parameters(2, 2)]);
% %         hold on;
% %         plot(x, Small);
%         if(value >= Domain(2, 1) && value <= Domain(2, 2))
%             membership_value(2) = Small(value - (Domain(2, 1)-1));
%         else
%             membership_value(2) = 0;
%         end
% 
%         x = Domain(3, 1):1:Domain(3, 2);
%         Medium = gaussmf(x,[Parameters(3, 1) Parameters(3, 2)]);
% %         hold on;
% %         plot(x, Medium);
%         if(value >= Domain(3, 1) && value <= Domain(3, 2))
%             membership_value(3) = Medium(value - (Domain(3, 1)-1));
%         else
%             membership_value(3) = 0;
%         end
% 
%         x = Domain(4, 1):1:Domain(4, 2);
%         Big = gaussmf(x,[Parameters(4, 1) Parameters(4, 2)]);
% %         hold on;
% %         plot(x, Big);
%         if(value >= Domain(4, 1) && value <= Domain(4, 2))
%             membership_value(4) = Big(value - (Domain(4, 1)-1));
%         else
%             membership_value(4) = 0;
%         end
% 
%         x = Domain(5, 1):1:Domain(5, 2);
%         Very_Big = smf(x,[Parameters(5, 1) Parameters(5, 2)]);
% %         hold on;
% %         plot(x, Very_Big);
%         if(value >= Domain(5, 1) && value <= Domain(5, 2))
%             membership_value(5) = Very_Big(value - (Domain(5, 1)-1));
%         else
%             membership_value(5) = 0;
%         end
        membership_value(1) = evalmf(value, [Parameters(1) Parameters(2)], 'zmf');
        membership_value(2) = evalmf(value, [Parameters(3) Parameters(4)], 'gaussmf');
        membership_value(3) = evalmf(value, [Parameters(5) Parameters(6)], 'gaussmf');
        membership_value(4) = evalmf(value, [Parameters(7) Parameters(8)], 'gaussmf');        
        membership_value(5) = evalmf(value, [Parameters(9) Parameters(10)], 'smf');
        return;
end



% if Eye Opening is Medium/Big              && Mouth Opening is Medium/Big/Very Big && Mouth Corner Displacement Positive Medium/Big          Then Happy  
% if Eye Opening is Medium/Big              && Mouth Opening is Very Small/Small && Mouth Corner Displacement Negative Small/Medium/Big       Then Sad
% if Eye Opening is Big                     && Mouth Opening is Very Small       && Mouth Corner Displacement Negative/Positive Small/Medium  Then Neutral
% if Eye Opening is Very Big                && Mouth Opening is Small/Medium/Big && Mouth Corner Displacement Negative/Positive Small/Medium  Then Fear
% if Eye Opening is Very Big                && Mouth Opening is Vrey Big         && Mouth Corner Displacement Positive Small                  Then Surprise
% if Eye Opening is Medium/Big              && Mouth Opening is Vrey Small       && Mouth Corner Displacement Negative Medium                 Then Angry
% if Eye Opening is Very Small/Small/Medium && Mouth Opening is very Small/Small && Mouth Corner Displacement Negative/Positive Small/Medium  Then Disgust


function [EmotionBelieve Emotion_Intencity Max_Believe] = Fuzzy_Emotion_Recognition(MouthOpening, MouthCornerDisplacement, EyeOpening, EyebrowConstriction, NoseSideWrinkle, MouthLength, MeanIntensity, EyebrowSlope, SadBelieve)

happy_believe1 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe2 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe3 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe4 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe5 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe6 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe7 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe8 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe9 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe10 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe11 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe12 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe13 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe14 = min([EyeOpening(2) MouthOpening(4) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe15 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe16 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe17 = min([EyeOpening(2) MouthOpening(4) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe18 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe19 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe20 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe21 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe22 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe23 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe24 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe25 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe26 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe27 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe28 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(3) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe29 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe30 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(2) min(1-EyebrowConstriction(1), max(EyebrowConstriction)) MouthLength(2)]);
happy_believe31 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(2) MouthLength(2)]);
happy_believe32 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(2) MouthLength(2)]);

very_happy = max([happy_believe5 happy_believe6 happy_believe7 happy_believe8 happy_believe9 happy_believe10 happy_believe11 happy_believe12 happy_believe13 happy_believe14 happy_believe15 happy_believe18 happy_believe31]);
moderately_happy = max([happy_believe1 happy_believe2 happy_believe3 happy_believe4 happy_believe16 happy_believe17 happy_believe23 happy_believe24 happy_believe25 happy_believe26 happy_believe27 happy_believe28 happy_believe32]);
not_so_happy = max([happy_believe19 happy_believe20 happy_believe21 happy_believe22 happy_believe29 happy_believe30]);

% happy_believe = 0.33*very_happy + 0.33*moderately_happy + 0.33*not_so_happy;
[happy_believe happy_intensity] = max([very_happy moderately_happy not_so_happy]);
% [happy_believe hi] = max([happy_believe1 happy_believe2 happy_believe3 happy_believe4 happy_believe5 happy_believe6 happy_believe7 happy_believe8 happy_believe9 happy_believe10 happy_believe11 happy_believe12 happy_believe13 happy_believe14 happy_believe15 happy_believe16 happy_believe17 happy_believe18 happy_believe19 happy_believe20 happy_believe21 happy_believe22 happy_believe23 happy_believe24 happy_believe25 happy_believe26 happy_believe27 happy_believe28 happy_believe29 happy_believe30 happy_believe31 happy_believe32]);


sad_believe1 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(3)]);%1
sad_believe2 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(3)]);
% sad_believe3 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(3)]);
% sad_believe4 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(3)]);
sad_believe5 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(3) EyebrowSlope(2)]);%2
sad_believe6 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(3) EyebrowSlope(2)]);
% sad_believe7 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(3)]);
% sad_believe8 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(3)]);
sad_believe9 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(3)]);%3%%
sad_believe10 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(3)]);%4%%
% sad_believe11 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(3)]);
% sad_believe12 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(3)]);
% sad_believe13 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(3)]);%5%%
sad_believe14 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(3) EyebrowSlope(2)]);%6%%
% sad_believe15 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(3)]);
% sad_believe16 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(3)]);
sad_believe17 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(2)]);%7
sad_believe18 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(2)]);
% sad_believe19 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(2)]);
% sad_believe20 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(2)]);
sad_believe21 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(2) EyebrowSlope(2)]);%8
sad_believe22 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(2)]);
% sad_believe23 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(2)]);
% sad_believe24 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(2)]);
sad_believe25 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(2)]);%9%%
sad_believe26 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(2)]);%10%%
% sad_believe27 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(2)]);
% sad_believe28 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(2)]);
sad_believe29 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(2) EyebrowSlope(2)]);%11%%
sad_believe30 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(2) EyebrowSlope(2)]);%12%%
% sad_believe31 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(2)]);
% sad_believe32 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(2)]);
sad_believe33 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
sad_believe34 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(1)]);
% sad_believe35 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
% sad_believe36 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(1)]);
sad_believe37 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
sad_believe38 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(1) EyebrowSlope(2)]);
% sad_believe39 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
% sad_believe40 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(1) EyebrowSlope(2)]);
sad_believe41 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(1)]);
sad_believe42 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(1)]);
% sad_believe43 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(1)]);
% sad_believe44 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(1)]);
sad_believe45 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(1) EyebrowSlope(2)]);
% % % % sad_believe46 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(1)]);
% sad_believe47 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(1)]);
% sad_believe48 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(1)]);
sad_believe49 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(3)]);
sad_believe50 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(3)]);
% sad_believe51 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(3)]);
% sad_believe52 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(3)]);
sad_believe53 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(3)]);
sad_believe54 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(3)]);
% sad_believe55 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(3)]);
% sad_believe56 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(3)]);
sad_believe57 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(2)]);
sad_believe58 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(2)]);
% sad_believe59 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(2)]);
% sad_believe60 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(2)]);
sad_believe61 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(2)]);
sad_believe62 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(2)]);
% sad_believe63 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(2)]);
% sad_believe64 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(2)]);
sad_believe65 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
sad_believe66 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(1)]);
% sad_believe67 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
% sad_believe68 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(1)]);
sad_believe69 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(1)]);
sad_believe70 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(1)]);
% sad_believe71 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(1)]);
% sad_believe72 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(1)]);
sad_believe73 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(1) EyebrowSlope(2)]);
sad_believe74 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(1) EyebrowSlope(2)]);
sad_believe75 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(1) EyebrowSlope(2)]);
sad_believe76 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(2) EyebrowSlope(2)]);
sad_believe77 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(2) EyebrowSlope(2)]);
sad_believe78 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(2) EyebrowSlope(2)]);
sad_believe79 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(3) EyebrowSlope(2)]);
sad_believe80 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(3) EyebrowSlope(2)]);
sad_believe81 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(3) EyebrowSlope(2)]);
sad_believe82 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(3) EyebrowSlope(2)]);
sad_believe83 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4) EyebrowSlope(2) SadBelieve(2)]);

very_sad = max([sad_believe1 sad_believe2 sad_believe5 sad_believe6 sad_believe17 sad_believe18 sad_believe21 sad_believe22 sad_believe29 sad_believe30 sad_believe33 sad_believe34 sad_believe37 sad_believe38 sad_believe45 sad_believe49 sad_believe50 sad_believe57 sad_believe58 sad_believe65 sad_believe66 sad_believe69 sad_believe70 sad_believe73 sad_believe74 sad_believe75 sad_believe76 sad_believe77 sad_believe78 sad_believe80 sad_believe81 sad_believe82 sad_believe83]); 
moderately_sad = max([sad_believe9 sad_believe10 sad_believe14 sad_believe25 sad_believe26 sad_believe41 sad_believe42 sad_believe61 sad_believe62 sad_believe79]);
not_so_sad = max([sad_believe53 sad_believe54]);

% sad_believe = 0.33*very_sad + 0.33*moderately_sad + 0.33*not_so_sad;
[sad_believe sad_intensity] = max([very_sad moderately_sad not_so_sad]);
% [sad_believe si] = max([sad_believe1 sad_believe2 sad_believe5 sad_believe6 sad_believe9 sad_believe10 sad_believe14 sad_believe17 sad_believe18 sad_believe21 sad_believe22 sad_believe25 sad_believe26 sad_believe29 sad_believe30 sad_believe33 sad_believe34 sad_believe37 sad_believe38 sad_believe41 sad_believe42 sad_believe45 sad_believe49 sad_believe50 sad_believe53 sad_believe54 sad_believe57 sad_believe58 sad_believe61 sad_believe62 sad_believe65 sad_believe66 sad_believe69 sad_believe70 sad_believe73 sad_believe74 sad_believe75 sad_believe76 sad_believe77 sad_believe78 sad_believe79 sad_believe80 sad_believe81 sad_believe82 sad_believe83]);




fear_believe1 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe2 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe3 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe4 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe5 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe6 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe7 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(1) EyebrowSlope(1)]);
fear_believe8 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(2) 1-NoseSideWrinkle EyebrowSlope(1)]);
fear_believe9 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(1)  EyebrowConstriction(3) 1-NoseSideWrinkle EyebrowSlope(1)]);
fear_believe10 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe11 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe12 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe13 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(1)]);
fear_believe14 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe15 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe16 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle]);
fear_believe17 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe18 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe19 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(1) EyebrowSlope(1)]);
fear_believe20 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle EyebrowSlope(1)]);
fear_believe21 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(3) 1-NoseSideWrinkle EyebrowSlope(1)]);
fear_believe22 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(3) 1-NoseSideWrinkle]);
fear_believe23 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe24 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe25 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
fear_believe26 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe27 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe28 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(2) 1-NoseSideWrinkle]);
fear_believe29 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(2) 1-NoseSideWrinkle]);
fear_believe30 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
% fear_believe31 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
% fear_believe32 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(2) 1-NoseSideWrinkle]);
% fear_believe33 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(3) 1-NoseSideWrinkle]);
fear_believe34 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(3) 1-NoseSideWrinkle]);
fear_believe35 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe36 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe37 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(1)]);
fear_believe38 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe39 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe40 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(2) 1-NoseSideWrinkle]);
fear_believe41 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe42 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
% fear_believe43 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(1)]);
% fear_believe44 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(2) 1-NoseSideWrinkle]);
% fear_believe45 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(3) 1-NoseSideWrinkle]);
fear_believe46 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(3) 1-NoseSideWrinkle]);
fear_believe47 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe48 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe49 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe50 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe51 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe52 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe53 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(3) 1-NoseSideWrinkle EyebrowSlope(1)]);
fear_believe54 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe55 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe56 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe57 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe58 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe59 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe60 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe61 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(1)]);
fear_believe62 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe63 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe64 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe65 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe66 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe67 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe68 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe69 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe70 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe71 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe72 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe73 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe74 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe75 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe76 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe77 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe78 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe79 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe80 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe81 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe82 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(1) MeanIntensity(1)]);
fear_believe83 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe84 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(3) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe85 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(1)]);
fear_believe86 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle]);
fear_believe87 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(2) EyebrowConstriction(1)]);
fear_believe88 = min([EyeOpening(5) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(3)]);
fear_believe89 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(1) 1-NoseSideWrinkle MeanIntensity(1) EyebrowSlope(1)]);
fear_believe90 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe91 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(1) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe92 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe93 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(2) MeanIntensity(1)]);
fear_believe94 = min([EyeOpening(5) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(1)]);
fear_believe95 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(2) EyebrowConstriction(3) MeanIntensity(1)]);

very_fear = max([fear_believe1 fear_believe2 fear_believe3 fear_believe4 fear_believe5 fear_believe6 fear_believe7 fear_believe8 fear_believe9 fear_believe10 fear_believe11 fear_believe12 fear_believe13 fear_believe14 fear_believe15 fear_believe16 fear_believe17 fear_believe18 fear_believe19 fear_believe20 fear_believe21 fear_believe22 fear_believe23 fear_believe24 fear_believe25 fear_believe26 fear_believe27 fear_believe28 fear_believe29 fear_believe30 fear_believe34 fear_believe35 fear_believe36 fear_believe37 fear_believe38 fear_believe39 fear_believe40 fear_believe41 fear_believe42 fear_believe46 fear_believe47 fear_believe48 fear_believe49 fear_believe50 fear_believe51 fear_believe52 fear_believe55 fear_believe57 fear_believe58 fear_believe59  fear_believe60 fear_believe61 fear_believe62 fear_believe63 fear_believe64 fear_believe65 fear_believe66 fear_believe67 fear_believe68 fear_believe69 fear_believe70 fear_believe71 fear_believe72 fear_believe73 fear_believe74 fear_believe75 fear_believe76 fear_believe77 fear_believe78 fear_believe79 fear_believe80 fear_believe81 fear_believe82 fear_believe83 fear_believe84 fear_believe85 fear_believe86 fear_believe87 fear_believe88 fear_believe89 fear_believe90 fear_believe91 fear_believe92 fear_believe93 fear_believe94 fear_believe95]);
moderately_fear = max([fear_believe54 fear_believe56]);
not_so_fear = fear_believe53;

% fear_believe = 0.33*very_fear + 0.33*moderately_fear + 0.33*not_so_fear;
[fear_believe fear_intensity] = max([very_fear moderately_fear not_so_fear]);
% fear_believe = max([fear_believe1 fear_believe2 fear_believe3 fear_believe4 fear_believe5 fear_believe6 fear_believe7 fear_believe8 fear_believe9 fear_believe10 fear_believe11 fear_believe12 fear_believe13 fear_believe14 fear_believe15 fear_believe16 fear_believe17 fear_believe18 fear_believe19 fear_believe20 fear_believe21 fear_believe22 fear_believe23 fear_believe24 fear_believe25 fear_believe26 fear_believe27 fear_believe28 fear_believe29 fear_believe30 fear_believe31 fear_believe32 fear_believe33 fear_believe34 fear_believe35 fear_believe36 fear_believe37 fear_believe38 fear_believe39 fear_believe40 fear_believe41 fear_believe42 fear_believe43 fear_believe44 fear_believe45 fear_believe46 fear_believe47 fear_believe48]);
% [fear_believe fi] = max([fear_believe1 fear_believe2 fear_believe3 fear_believe4 fear_believe5 fear_believe6 fear_believe7 fear_believe8 fear_believe9 fear_believe10 fear_believe11 fear_believe12 fear_believe13 fear_believe14 fear_believe15 fear_believe16 fear_believe17 fear_believe18 fear_believe19 fear_believe20 fear_believe21 fear_believe22 fear_believe23 fear_believe24 fear_believe25 fear_believe26 fear_believe27 fear_believe28 fear_believe29 fear_believe30 fear_believe34 fear_believe35 fear_believe36 fear_believe37 fear_believe38 fear_believe39 fear_believe40 fear_believe41 fear_believe42 fear_believe46 fear_believe47 fear_believe48 fear_believe49 fear_believe50 fear_believe51 fear_believe52 fear_believe53 fear_believe54 fear_believe55 fear_believe56 fear_believe57 fear_believe58 fear_believe59 fear_believe60 fear_believe61 fear_believe62 fear_believe63 fear_believe64 fear_believe65 fear_believe66 fear_believe67 fear_believe68 fear_believe69 fear_believe70 fear_believe71 fear_believe72 fear_believe73 fear_believe74 fear_believe75 fear_believe76 fear_believe77 fear_believe78 fear_believe79 fear_believe80 fear_believe81 fear_believe82 fear_believe83 fear_believe84 fear_believe85 fear_believe86 fear_believe87 fear_believe88 fear_believe89 fear_believe90 fear_believe91 fear_believe92 fear_believe93 fear_believe94 fear_believe95]);




surprise_believe1 = min([EyeOpening(5) MouthOpening(5) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1)]);
surprise_believe2 = min([EyeOpening(5) MouthOpening(5) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe3 = min([EyeOpening(5) MouthOpening(5) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1)]);
surprise_believe4 = min([EyeOpening(5) MouthOpening(5) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe5 = min([EyeOpening(4) MouthOpening(5) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1)]);
surprise_believe6 = min([EyeOpening(4) MouthOpening(5) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe7 = min([EyeOpening(4) MouthOpening(5) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1)]);
surprise_believe8 = min([EyeOpening(4) MouthOpening(5) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe9 = min([EyeOpening(5) MouthOpening(4) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe10 = min([EyeOpening(5) MouthOpening(4) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe11 = min([EyeOpening(5) MouthOpening(4) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe12 = min([EyeOpening(5) MouthOpening(4) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe13 = min([EyeOpening(4) MouthOpening(4) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe14 = min([EyeOpening(4) MouthOpening(4) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe15 = min([EyeOpening(4) MouthOpening(4) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe16 = min([EyeOpening(4) MouthOpening(4) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe17 = min([EyeOpening(5) MouthOpening(2) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe18 = min([EyeOpening(4) MouthOpening(2) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(1) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe19 = min([EyeOpening(4) MouthOpening(2) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe20 = min([EyeOpening(5) MouthOpening(2) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(3) 1-NoseSideWrinkle MouthLength(1) MeanIntensity(2)]);
surprise_believe21 = min([EyeOpening(5) MouthOpening(3) min([1-MouthCornerDisplacement(6) max(MouthCornerDisplacement)]) EyebrowConstriction(2) 1-NoseSideWrinkle MeanIntensity(2)]);

very_surprise = max([surprise_believe1 surprise_believe2 surprise_believe3 surprise_believe4 surprise_believe5 surprise_believe6 surprise_believe7 surprise_believe8  surprise_believe9 surprise_believe10 surprise_believe11 surprise_believe12 surprise_believe13 surprise_believe14 surprise_believe15 surprise_believe16 surprise_believe21]);
moderately_surprise = max([surprise_believe17 surprise_believe18 surprise_believe19 surprise_believe20]);

% surprise_believe = 0.33*very_surprise + 0.33*moderately_surprise;
[surprise_believe surprise_intensity] = max([very_surprise moderately_surprise]);
% [surprise_believe sui] = max([surprise_believe1 surprise_believe2 surprise_believe3 surprise_believe4 surprise_believe5 surprise_believe6 surprise_believe7 surprise_believe8 surprise_believe9 surprise_believe10 surprise_believe11 surprise_believe12 surprise_believe13 surprise_believe14 surprise_believe15 surprise_believe16 surprise_believe17 surprise_believe18 surprise_believe19 surprise_believe20 surprise_believe21]);



angry_believe1 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
angry_believe2 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
angry_believe3 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
angry_believe4 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
angry_believe5 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
angry_believe6 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
angry_believe7 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
angry_believe8 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4) SadBelieve(1)]);
angry_believe9 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
angry_believe10 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4) SadBelieve(1)]);
angry_believe11 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(4)]);
angry_believe12 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(4) SadBelieve(1)]);
angry_believe13 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
angry_believe14 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
angry_believe15 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
angry_believe16 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(4) SadBelieve(1)]);
angry_believe17 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
angry_believe18 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4) SadBelieve(1)]);
angry_believe19 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
angry_believe20 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4) SadBelieve(1)]);
angry_believe21 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
angry_believe22 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(4) SadBelieve(1)]);
angry_believe23 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
angry_believe24 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(4) SadBelieve(1)]);
% angry_believe25 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
% angry_believe26 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
% angry_believe27 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
% angry_believe28 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
% angry_believe29 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
% angry_believe30 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
% angry_believe31 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
% angry_believe32 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
% angry_believe33 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
% angry_believe34 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
% angry_believe35 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(4)]);
% angry_believe36 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(4)]);
% angry_believe37 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
% angry_believe38 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
% angry_believe39 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
% angry_believe40 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
% angry_believe41 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
% angry_believe42 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
% angry_believe43 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
% angry_believe44 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
% angry_believe45 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
% angry_believe46 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(4)]);
% angry_believe47 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
% angry_believe48 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
% angry_believe49 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
% angry_believe50 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
% angry_believe51 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
% angry_believe52 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
% angry_believe53 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
% angry_believe54 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
% angry_believe55 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
% angry_believe56 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
% angry_believe57 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
% angry_believe58 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
% angry_believe59 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
% angry_believe60 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
% angry_believe61 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
angry_believe62 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4) SadBelieve(1)]);

very_angry = max([angry_believe1 angry_believe2 angry_believe3 angry_believe4 angry_believe5 angry_believe6 angry_believe7 angry_believe8 angry_believe9 angry_believe11 angry_believe13 angry_believe14 angry_believe15 angry_believe16 angry_believe17 angry_believe18 angry_believe19 angry_believe21 angry_believe23 angry_believe62]);
moderately_angry = max([angry_believe10 angry_believe12 angry_believe20 angry_believe22 angry_believe24]);

% angry_believe = 0.33*very_angry + 0.33*moderately_angry;
[angry_believe angry_intensity] = max([very_angry moderately_angry]);
% [angry_believe ai] = max([angry_believe1 angry_believe2 angry_believe3 angry_believe4 angry_believe5 angry_believe6 angry_believe7 angry_believe8 angry_believe9 angry_believe10 angry_believe11 angry_believe12 angry_believe13 angry_believe14 angry_believe15 angry_believe16 angry_believe17 angry_believe18 angry_believe19 angry_believe20 angry_believe21 angry_believe22 angry_believe23 angry_believe24 angry_believe25 angry_believe26 angry_believe27 angry_believe28 angry_believe29 angry_believe30 angry_believe31 angry_believe32 angry_believe33 angry_believe34 angry_believe35 angry_believe36 angry_believe37 angry_believe38 angry_believe39 angry_believe40 angry_believe41 angry_believe42 angry_believe43 angry_believe44 angry_believe45 angry_believe46 angry_believe47 angry_believe48 angry_believe49 angry_believe50 angry_believe51 angry_believe52 angry_believe53 angry_believe54 angry_believe55 angry_believe56 angry_believe57 angry_believe58 angry_believe59 angry_believe60 angry_believe61]);
% [angry_believe ai] = max([angry_believe1 angry_believe2 angry_believe3 angry_believe4 angry_believe5 angry_believe6 angry_believe7 angry_believe8 angry_believe9 angry_believe10 angry_believe11 angry_believe12 angry_believe13 angry_believe14 angry_believe15 angry_believe16 angry_believe17 angry_believe18 angry_believe19 angry_believe20 angry_believe21 angry_believe22 angry_believe23 angry_believe24 angry_believe62]);



disgust_believe1 = min([EyeOpening(1) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe2 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe3 = min([EyeOpening(1) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe4 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe5 = min([EyeOpening(1) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe6 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe7 = min([EyeOpening(1) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe8 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe9 = min([EyeOpening(1) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe10 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe11 = min([EyeOpening(1) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe12 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe13 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe14 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe15 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe16 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe17 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe18 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe19 = min([EyeOpening(1) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe20 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe21 = min([EyeOpening(1) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe22 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe23 = min([EyeOpening(1) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe24 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe25 = min([EyeOpening(1) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe26 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe27 = min([EyeOpening(1) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe28 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe29 = min([EyeOpening(1) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe30 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe31 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe32 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe33 = min([EyeOpening(3) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe34 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe35 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe36 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe37 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe38 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe39 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe40 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe41 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe42 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe43 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe44 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe45 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe46 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe47 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe48 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe49 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe50 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe51 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe52 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe53 = min([EyeOpening(1) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
disgust_believe54 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
disgust_believe55 = min([EyeOpening(1) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
disgust_believe56 = min([EyeOpening(2) MouthOpening(2) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
disgust_believe57 = min([EyeOpening(1) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
disgust_believe58 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
disgust_believe59 = min([EyeOpening(1) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
disgust_believe60 = min([EyeOpening(2) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
disgust_believe61 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
disgust_believe62 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
disgust_believe63 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(4)]);
disgust_believe64 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(6) EyebrowConstriction(5)]);
disgust_believe65 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe66 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe67 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe68 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe69 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe70 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe71 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe72 = min([EyeOpening(4) MouthOpening(2) MouthCornerDisplacement(5) EyebrowConstriction(5)]);
disgust_believe73 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe74 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe75 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(5) EyebrowConstriction(4)]);
disgust_believe76 = min([EyeOpening(4) MouthOpening(3) MouthCornerDisplacement(1) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe77 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe78 = min([EyeOpening(1) MouthOpening(3) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe79 = min([EyeOpening(1) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe80 = min([EyeOpening(1) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe81 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe82 = min([EyeOpening(2) MouthOpening(3) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe83 = min([EyeOpening(2) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe84 = min([EyeOpening(2) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe85 = min([EyeOpening(1) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe86 = min([EyeOpening(1) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe87 = min([EyeOpening(2) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe88 = min([EyeOpening(2) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe89 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe90 = min([EyeOpening(3) MouthOpening(3) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe91 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe92 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe93 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe94 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe95 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe96 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(5)]);
disgust_believe97 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe98 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe99 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4) NoseSideWrinkle MouthLength(1)]);
disgust_believe100 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe101 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle MouthLength(1)]);
disgust_believe102 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe103 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe104 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe105 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe106 = min([EyeOpening(3) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe107 = min([EyeOpening(4) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(4) NoseSideWrinkle MouthLength(1)]);
disgust_believe108 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe109 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe110 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe111 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe112 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe113 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe114 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe115 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe116 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe117 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe118 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe119 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe120 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe121 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe122 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe123 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe124 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe125 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(5) EyebrowConstriction(4) NoseSideWrinkle MouthLength(1)]);
disgust_believe126 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe127 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe128 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe129 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(1) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe130 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe131 = min([EyeOpening(5) MouthOpening(1) MouthCornerDisplacement(6) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe132 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe133 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe134 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe135 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe136 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe137 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe138 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe139 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe140 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe141 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe142 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe143 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe144 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe145 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe146 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(5) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe147 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe148 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe149 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe150 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(4) NoseSideWrinkle]);
disgust_believe151 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(6) EyebrowConstriction(5) NoseSideWrinkle]);
disgust_believe152 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(4) NoseSideWrinkle MouthLength(1)]);
disgust_believe153 = min([EyeOpening(1) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(5) NoseSideWrinkle MouthLength(1)]);
disgust_believe154 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(4) NoseSideWrinkle MouthLength(1)]);
disgust_believe155 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(5) NoseSideWrinkle MouthLength(1)]);
disgust_believe156 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(4) NoseSideWrinkle MouthLength(1)]);
disgust_believe157 = min([EyeOpening(3) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(5) NoseSideWrinkle MouthLength(1)]);
disgust_believe158 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(4)]);
disgust_believe159 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(5)]);
disgust_believe160 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe161 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe162 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe163 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(3) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe164 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(3) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe165 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe166 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe167 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(3) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe168 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(3) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe169 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe170 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe171 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe172 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe173 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe174 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(3) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe175 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(3) EyebrowConstriction(4) MouthLength(1)]);
disgust_believe176 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe177 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(2) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe178 = min([EyeOpening(5) MouthOpening(4) MouthCornerDisplacement(3) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe179 = min([EyeOpening(5) MouthOpening(5) MouthCornerDisplacement(3) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe180 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(5) MouthLength(1)]);
disgust_believe181 = min([EyeOpening(2) MouthOpening(5) MouthCornerDisplacement(1) EyebrowConstriction(2) MouthLength(1)]);
disgust_believe182 = min([EyeOpening(3) MouthOpening(4) MouthCornerDisplacement(4) EyebrowConstriction(2) MouthLength(1)]);
disgust_believe183 = min([EyeOpening(4) MouthOpening(4) MouthCornerDisplacement(1) EyebrowConstriction(2) NoseSideWrinkle MouthLength(1)]);
disgust_believe184 = min([EyeOpening(4) MouthOpening(5) MouthCornerDisplacement(4) EyebrowConstriction(3) NoseSideWrinkle MouthLength(1)]);

very_disgust = max([disgust_believe1 disgust_believe2 disgust_believe3 disgust_believe4 disgust_believe5 disgust_believe6 disgust_believe7 disgust_believe8 disgust_believe9 disgust_believe10 disgust_believe11 disgust_believe12 disgust_believe13 disgust_believe14 disgust_believe15 disgust_believe16 disgust_believe17 disgust_believe18 disgust_believe19 disgust_believe20 disgust_believe21 disgust_believe22 disgust_believe23 disgust_believe24 disgust_believe25 disgust_believe26 disgust_believe27 disgust_believe28 disgust_believe29 disgust_believe30 disgust_believe31 disgust_believe32 disgust_believe33 disgust_believe34 disgust_believe35 disgust_believe36 disgust_believe37 disgust_believe38 disgust_believe39 disgust_believe40 disgust_believe41 disgust_believe42 disgust_believe43 disgust_believe46 disgust_believe47 disgust_believe48 disgust_believe49 disgust_believe50 disgust_believe51 disgust_believe52 disgust_believe53 disgust_believe54 disgust_believe54 disgust_believe55 disgust_believe56 disgust_believe57 disgust_believe58 disgust_believe58 disgust_believe59 disgust_believe60 disgust_believe61 disgust_believe62 disgust_believe63 disgust_believe64 disgust_believe68 disgust_believe69 disgust_believe70 disgust_believe71 disgust_believe72 disgust_believe73 disgust_believe74 disgust_believe75 disgust_believe76 disgust_believe77 disgust_believe78 disgust_believe79 disgust_believe80 disgust_believe81 disgust_believe82 disgust_believe83 disgust_believe84 disgust_believe85 disgust_believe86 disgust_believe87 disgust_believe88 disgust_believe89 disgust_believe90 disgust_believe91 disgust_believe92 disgust_believe93 disgust_believe94 disgust_believe95 disgust_believe96 disgust_believe97 disgust_believe98 disgust_believe99 disgust_believe100 disgust_believe101 disgust_believe102 disgust_believe103 disgust_believe104 disgust_believe105 disgust_believe106 disgust_believe107 disgust_believe108 disgust_believe109 disgust_believe110 disgust_believe111 disgust_believe112 disgust_believe113 disgust_believe114 disgust_believe115 disgust_believe116 disgust_believe117 disgust_believe118 disgust_believe119 disgust_believe120 disgust_believe121 disgust_believe122 disgust_believe123 disgust_believe124 disgust_believe125 disgust_believe126 disgust_believe127 disgust_believe128 disgust_believe129 disgust_believe130 disgust_believe131 disgust_believe132 disgust_believe133 disgust_believe134 disgust_believe135 disgust_believe136 disgust_believe137 disgust_believe138 disgust_believe139 disgust_believe140 disgust_believe141 disgust_believe142 disgust_believe143 disgust_believe144 disgust_believe145 disgust_believe146 disgust_believe147 disgust_believe148 disgust_believe149 disgust_believe150 disgust_believe151 disgust_believe152 disgust_believe153 disgust_believe154 disgust_believe155 disgust_believe156 disgust_believe157 disgust_believe158 disgust_believe159 disgust_believe160 disgust_believe161 disgust_believe162 disgust_believe163 disgust_believe164 disgust_believe165 disgust_believe166 disgust_believe167 disgust_believe168 disgust_believe169 disgust_believe170 disgust_believe171 disgust_believe172 disgust_believe172 disgust_believe173 disgust_believe174 disgust_believe175 disgust_believe176 disgust_believe177 disgust_believe178 disgust_believe179 disgust_believe180 disgust_believe181 disgust_believe182 disgust_believe183  disgust_believe184]);
moderately_disgust = max([disgust_believe44 disgust_believe45 disgust_believe65 disgust_believe67]);
not_so_disgust = disgust_believe66;

% disgust_believe = 0.33*very_disgust + 0.33*moderately_disgust + 0.33*not_so_disgust;
[disgust_believe disgust_intensity] = max([very_disgust moderately_disgust not_so_disgust]);
% [disgust_believe di] = max([disgust_believe1 disgust_believe2 disgust_believe3 disgust_believe4 disgust_believe5 disgust_believe6 disgust_believe7 disgust_believe8 disgust_believe9 disgust_believe10 disgust_believe11 disgust_believe12 disgust_believe13 disgust_believe14 disgust_believe15 disgust_believe16 disgust_believe17 disgust_believe18 disgust_believe19 disgust_believe20 disgust_believe21 disgust_believe22 disgust_believe23 disgust_believe24 disgust_believe25 disgust_believe26 disgust_believe27 disgust_believe28 disgust_believe29 disgust_believe30 disgust_believe31 disgust_believe32 disgust_believe33 disgust_believe34 disgust_believe35 disgust_believe36 disgust_believe37 disgust_believe38 disgust_believe39 disgust_believe40 disgust_believe41 disgust_believe42 disgust_believe43 disgust_believe44 disgust_believe45 disgust_believe46 disgust_believe47 disgust_believe48 disgust_believe49 disgust_believe50 disgust_believe51 disgust_believe52 disgust_believe53 disgust_believe54 disgust_believe54 disgust_believe55 disgust_believe56 disgust_believe57 disgust_believe58 disgust_believe58 disgust_believe59 disgust_believe60 disgust_believe61 disgust_believe62 disgust_believe63 disgust_believe64 disgust_believe65 disgust_believe66 disgust_believe67 disgust_believe68 disgust_believe69 disgust_believe70 disgust_believe71 disgust_believe72 disgust_believe73 disgust_believe74 disgust_believe75 disgust_believe76 disgust_believe77 disgust_believe78 disgust_believe79 disgust_believe80 disgust_believe81 disgust_believe82 disgust_believe83 disgust_believe84 disgust_believe85 disgust_believe86 disgust_believe87 disgust_believe88 disgust_believe89 disgust_believe90 disgust_believe91 disgust_believe92 disgust_believe93 disgust_believe94 disgust_believe95 disgust_believe96 disgust_believe97 disgust_believe98 disgust_believe99 disgust_believe100 disgust_believe101 disgust_believe102 disgust_believe103 disgust_believe104 disgust_believe105 disgust_believe106 disgust_believe107 disgust_believe108 disgust_believe109 disgust_believe110 disgust_believe111 disgust_believe112 disgust_believe113 disgust_believe114 disgust_believe115 disgust_believe116 disgust_believe117 disgust_believe118 disgust_believe119 disgust_believe120 disgust_believe121 disgust_believe122 disgust_believe123 disgust_believe124 disgust_believe125 disgust_believe126 disgust_believe127 disgust_believe128 disgust_believe129 disgust_believe130 disgust_believe131 disgust_believe132 disgust_believe133 disgust_believe134 disgust_believe135 disgust_believe136 disgust_believe137 disgust_believe138 disgust_believe139 disgust_believe140 disgust_believe141 disgust_believe142 disgust_believe143 disgust_believe144 disgust_believe145 disgust_believe146 disgust_believe147 disgust_believe148 disgust_believe149 disgust_believe150 disgust_believe151 disgust_believe152 disgust_believe153 disgust_believe154 disgust_believe155 disgust_believe156 disgust_believe157 disgust_believe158 disgust_believe159 disgust_believe160 disgust_believe161 disgust_believe162 disgust_believe163 disgust_believe164 disgust_believe165 disgust_believe166 disgust_believe167 disgust_believe168 disgust_believe169 disgust_believe170 disgust_believe171 disgust_believe172 disgust_believe172 disgust_believe173 disgust_believe174 disgust_believe175 disgust_believe176 disgust_believe177 disgust_believe178 disgust_believe179 disgust_believe180 disgust_believe181 disgust_believe182 disgust_believe183  disgust_believe184]);



EmotionBelieve = [happy_believe sad_believe fear_believe surprise_believe angry_believe disgust_believe];
[Max_Believe max_index] = max(EmotionBelieve);

if(length(Max_Believe) > 1)
    Emotion_Intencity = 'Blend Emotions';
else
    switch(max_index)
        case 1
            if(happy_intensity == 1)
                Emotion_Intencity = 'Very Happy';
            else
                if(happy_intensity == 2)
                    Emotion_Intencity = 'Moderately Happy';
                else
                    Emotion_Intencity = 'Not-So Happy';
                end
            end
            
        case 2
            if(sad_intensity == 1)
                Emotion_Intencity = 'Very Sad';
            else
                if(sad_intensity == 2)
                    Emotion_Intencity = 'Moderately Sad';
                else
                    Emotion_Intencity = 'Not-So Sad';
                end
            end

        case 3
            if(fear_intensity == 1)
                Emotion_Intencity = 'Very Fear';
            else
                if(fear_intensity == 2)
                    Emotion_Intencity = 'Moderately Fear';
                else
                    Emotion_Intencity = 'Not-So Fear';
                end
            end
            
        case 4
            if(surprise_intensity == 1)
                Emotion_Intencity = 'Very Surprise';
            else
                Emotion_Intencity = 'Moderately Surprise';
            end
            
        case 5
            if(angry_intensity == 1)
                Emotion_Intencity = 'Very Angry';
            else
                Emotion_Intencity = 'Moderately Angry';
            end
            
        case 6
            if(disgust_intensity == 1)
                Emotion_Intencity = 'Very Disgust';
            else
                if(disgust_intensity == 2)
                    Emotion_Intencity = 'Moderately Disgust';
                else
                    Emotion_Intencity = 'Not-So Disgust';
                end
            end

    end
end


function [EmotionBelieve Emotion_Intencity Max_Believe] = Fuzzy_Emotion_Recognition2(MouthOpening, MouthCornerDisplacement, EyeOpening, EyebrowConstriction, NoseSideWrinkle, MouthLength, MeanIntensity, EyebrowSlope, SadBelieve)

happy_believe1 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe2 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe3 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe4 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe5 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe6 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe7 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe8 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe9 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe10 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe11 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe12 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe13 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe14 = EyeOpening(2) * MouthOpening(4) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe15 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe16 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe17 = EyeOpening(2) * MouthOpening(4) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe18 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe19 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe20 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe21 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe22 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe23 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe24 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe25 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe26 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe27 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe28 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(3) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe29 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe30 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(2) * (1-EyebrowConstriction(1)) * max(EyebrowConstriction) * MouthLength(2);
happy_believe31 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(2) * MouthLength(2);
happy_believe32 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(2) * MouthLength(2);

very_happy = max([happy_believe5 happy_believe6 happy_believe7 happy_believe8 happy_believe9 happy_believe10 happy_believe11 happy_believe12 happy_believe13 happy_believe14 happy_believe15 happy_believe18 happy_believe31]);
moderately_happy = max([happy_believe1 happy_believe2 happy_believe3 happy_believe4 happy_believe16 happy_believe17 happy_believe23 happy_believe24 happy_believe25 happy_believe26 happy_believe27 happy_believe28 happy_believe32]);
not_so_happy = max([happy_believe19 happy_believe20 happy_believe21 happy_believe22 happy_believe29 happy_believe30]);

[happy_believe happy_intensity] = max([very_happy moderately_happy not_so_happy]);
% [happy_believe hi] = max([happy_believe1 happy_believe2 happy_believe3 happy_believe4 happy_believe5 happy_believe6 happy_believe7 happy_believe8 happy_believe9 happy_believe10 happy_believe11 happy_believe12 happy_believe13 happy_believe14 happy_believe15 happy_believe16 happy_believe17 happy_believe18 happy_believe19 happy_believe20 happy_believe21 happy_believe22 happy_believe23 happy_believe24 happy_believe25 happy_believe26 happy_believe27 happy_believe28 happy_believe29 happy_believe30 happy_believe31 happy_believe32;


sad_believe1 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(3);%1
sad_believe2 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(3);
% sad_believe3 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(3);
% sad_believe4 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(3);
sad_believe5 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * EyebrowSlope(2);%2
sad_believe6 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * EyebrowSlope(2);
% sad_believe7 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(3);
% sad_believe8 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(3);
sad_believe9 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(3);%3%%
sad_believe10 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(3);%4%%
% sad_believe11 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(3);
% sad_believe12 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(3);
% sad_believe13 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(3);%5%%
sad_believe14 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(3) * EyebrowSlope(2);%6%%
% sad_believe15 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(3);
% sad_believe16 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(3);
sad_believe17 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(2);%7
sad_believe18 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(2);
% sad_believe19 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(2);
% sad_believe20 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(2);
sad_believe21 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * EyebrowSlope(2);%8
sad_believe22 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(2);
% sad_believe23 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(2);
% sad_believe24 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(2);
sad_believe25 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(2);%9%%
sad_believe26 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(2);%10%%
% sad_believe27 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(2);
% sad_believe28 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(2);
sad_believe29 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * EyebrowSlope(2);%11%%
sad_believe30 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * EyebrowSlope(2);%12%%
% sad_believe31 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(2);
% sad_believe32 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(2);
sad_believe33 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
sad_believe34 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(1);
% sad_believe35 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
% sad_believe36 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(1);
sad_believe37 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
sad_believe38 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * EyebrowSlope(2);
% sad_believe39 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
% sad_believe40 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * EyebrowSlope(2);
sad_believe41 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(1);
sad_believe42 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(1);
% sad_believe43 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(1);
% sad_believe44 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(1);
sad_believe45 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(1) * EyebrowSlope(2);
% % % % sad_believe46 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(1);
% sad_believe47 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(1);
% sad_believe48 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(1);
sad_believe49 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(3);
sad_believe50 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(3);
% sad_believe51 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(3);
% sad_believe52 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(3);
sad_believe53 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(3);
sad_believe54 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(3);
% sad_believe55 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(3);
% sad_believe56 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(3);
sad_believe57 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(2);
sad_believe58 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(2);
% sad_believe59 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(2);
% sad_believe60 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(2);
sad_believe61 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(2);
sad_believe62 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(2);
% sad_believe63 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(2);
% sad_believe64 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(2);
sad_believe65 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
sad_believe66 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(1);
% sad_believe67 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
% sad_believe68 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(1);
sad_believe69 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(1);
sad_believe70 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(1);
% sad_believe71 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(1);
% sad_believe72 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(1);
sad_believe73 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(1) * EyebrowSlope(2);
sad_believe74 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(1) * EyebrowSlope(2);
sad_believe75 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * EyebrowSlope(2);
sad_believe76 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * EyebrowSlope(2);
sad_believe77 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * EyebrowSlope(2);
sad_believe78 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * EyebrowSlope(2);
sad_believe79 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(3) * EyebrowSlope(2);
sad_believe80 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * EyebrowSlope(2);
sad_believe81 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * EyebrowSlope(2);
sad_believe82 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * EyebrowSlope(2);
sad_believe83 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * EyebrowSlope(2) * SadBelieve(2);

very_sad = max([sad_believe1 sad_believe2 sad_believe5 sad_believe6 sad_believe17 sad_believe18 sad_believe21 sad_believe22 sad_believe29 sad_believe30 sad_believe33 sad_believe34 sad_believe37 sad_believe38 sad_believe45 sad_believe49 sad_believe50 sad_believe57 sad_believe58 sad_believe65 sad_believe66 sad_believe69 sad_believe70 sad_believe73 sad_believe74 sad_believe75 sad_believe76 sad_believe77 sad_believe78 sad_believe80 sad_believe81 sad_believe82 sad_believe83]); 
moderately_sad = max([sad_believe9 sad_believe10 sad_believe14 sad_believe25 sad_believe26 sad_believe41 sad_believe42 sad_believe61 sad_believe62 sad_believe79]);
not_so_sad = max([sad_believe53 sad_believe54]);

[sad_believe sad_intensity] = max([very_sad moderately_sad not_so_sad]);
% [sad_believe si] = max([sad_believe1 sad_believe2 sad_believe5 sad_believe6 sad_believe9 sad_believe10 sad_believe14 sad_believe17 sad_believe18 sad_believe21 sad_believe22 sad_believe25 sad_believe26 sad_believe29 sad_believe30 sad_believe33 sad_believe34 sad_believe37 sad_believe38 sad_believe41 sad_believe42 sad_believe45 sad_believe49 sad_believe50 sad_believe53 sad_believe54 sad_believe57 sad_believe58 sad_believe61 sad_believe62 sad_believe65 sad_believe66 sad_believe69 sad_believe70 sad_believe73 sad_believe74 sad_believe75 sad_believe76 sad_believe77 sad_believe78 sad_believe79 sad_believe80 sad_believe81 sad_believe82 sad_believe83;



fear_believe1 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe2 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe3 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe4 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe5 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe6 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe7 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(1) * EyebrowSlope(1);
fear_believe8 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * EyebrowSlope(1);
fear_believe9 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(1) *  EyebrowConstriction(3) * (1-NoseSideWrinkle) * EyebrowSlope(1);
fear_believe10 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe11 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe12 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe13 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(1);
fear_believe14 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe15 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe16 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle);
fear_believe17 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe18 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe19 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(1) * EyebrowSlope(1);
fear_believe20 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * EyebrowSlope(1);
fear_believe21 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * EyebrowSlope(1);
fear_believe22 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(3) * (1-NoseSideWrinkle);
fear_believe23 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe24 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe25 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
fear_believe26 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe27 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe28 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * (1-NoseSideWrinkle);
fear_believe29 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * (1-NoseSideWrinkle);
fear_believe30 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
% fear_believe31 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
% fear_believe32 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * (1-NoseSideWrinkle);
% fear_believe33 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * (1-NoseSideWrinkle);
fear_believe34 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * (1-NoseSideWrinkle);
fear_believe35 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe36 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe37 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(1);
fear_believe38 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe39 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe40 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * (1-NoseSideWrinkle);
fear_believe41 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe42 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
% fear_believe43 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(1);
% fear_believe44 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * (1-NoseSideWrinkle);
% fear_believe45 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * (1-NoseSideWrinkle);
fear_believe46 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * (1-NoseSideWrinkle);
fear_believe47 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe48 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe49 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe50 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe51 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe52 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe53 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * EyebrowSlope(1);
fear_believe54 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe55 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe56 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe57 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe58 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe59 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe60 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe61 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(1);
fear_believe62 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe63 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe64 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe65 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe66 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe67 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe68 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe69 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe70 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe71 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe72 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe73 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe74 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe75 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe76 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe77 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe78 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe79 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe80 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe81 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe82 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(1) * MeanIntensity(1);
fear_believe83 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe84 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(3) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe85 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(1);
fear_believe86 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle);
fear_believe87 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(2) * EyebrowConstriction(1);
fear_believe88 = EyeOpening(5) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(3);
fear_believe89 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(1) * (1-NoseSideWrinkle) * MeanIntensity(1) * EyebrowSlope(1);
fear_believe90 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe91 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(1) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe92 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe93 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * MeanIntensity(1);
fear_believe94 = EyeOpening(5) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * (1-NoseSideWrinkle) * MeanIntensity(1);
fear_believe95 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(2) * EyebrowConstriction(3) * MeanIntensity(1);

very_fear = max([fear_believe1 fear_believe2 fear_believe3 fear_believe4 fear_believe5 fear_believe6 fear_believe7 fear_believe8 fear_believe9 fear_believe10 fear_believe11 fear_believe12 fear_believe13 fear_believe14 fear_believe15 fear_believe16 fear_believe17 fear_believe18 fear_believe19 fear_believe20 fear_believe21 fear_believe22 fear_believe23 fear_believe24 fear_believe25 fear_believe26 fear_believe27 fear_believe28 fear_believe29 fear_believe30 fear_believe34 fear_believe35 fear_believe36 fear_believe37 fear_believe38 fear_believe39 fear_believe40 fear_believe41 fear_believe42 fear_believe46 fear_believe47 fear_believe48 fear_believe49 fear_believe50 fear_believe51 fear_believe52 fear_believe55 fear_believe57 fear_believe58 fear_believe59  fear_believe60 fear_believe61 fear_believe62 fear_believe63 fear_believe64 fear_believe65 fear_believe66 fear_believe67 fear_believe68 fear_believe69 fear_believe70 fear_believe71 fear_believe72 fear_believe73 fear_believe74 fear_believe75 fear_believe76 fear_believe77 fear_believe78 fear_believe79 fear_believe80 fear_believe81 fear_believe82 fear_believe83 fear_believe84 fear_believe85 fear_believe86 fear_believe87 fear_believe88 fear_believe89 fear_believe90 fear_believe91 fear_believe92 fear_believe93 fear_believe94 fear_believe95]);
moderately_fear = max([fear_believe54 fear_believe56]);
not_so_fear = fear_believe53;

[fear_believe fear_intensity] = max([very_fear moderately_fear not_so_fear]);
% fear_believe = max([fear_believe1 fear_believe2 fear_believe3 fear_believe4 fear_believe5 fear_believe6 fear_believe7 fear_believe8 fear_believe9 fear_believe10 fear_believe11 fear_believe12 fear_believe13 fear_believe14 fear_believe15 fear_believe16 fear_believe17 fear_believe18 fear_believe19 fear_believe20 fear_believe21 fear_believe22 fear_believe23 fear_believe24 fear_believe25 fear_believe26 fear_believe27 fear_believe28 fear_believe29 fear_believe30 fear_believe31 fear_believe32 fear_believe33 fear_believe34 fear_believe35 fear_believe36 fear_believe37 fear_believe38 fear_believe39 fear_believe40 fear_believe41 fear_believe42 fear_believe43 fear_believe44 fear_believe45 fear_believe46 fear_believe47 fear_believe48;
% [fear_believe fi] = max([fear_believe1 fear_believe2 fear_believe3 fear_believe4 fear_believe5 fear_believe6 fear_believe7 fear_believe8 fear_believe9 fear_believe10 fear_believe11 fear_believe12 fear_believe13 fear_believe14 fear_believe15 fear_believe16 fear_believe17 fear_believe18 fear_believe19 fear_believe20 fear_believe21 fear_believe22 fear_believe23 fear_believe24 fear_believe25 fear_believe26 fear_believe27 fear_believe28 fear_believe29 fear_believe30 fear_believe34 fear_believe35 fear_believe36 fear_believe37 fear_believe38 fear_believe39 fear_believe40 fear_believe41 fear_believe42 fear_believe46 fear_believe47 fear_believe48 fear_believe49 fear_believe50 fear_believe51 fear_believe52 fear_believe53 fear_believe54 fear_believe55 fear_believe56 fear_believe57 fear_believe58 fear_believe59 fear_believe60 fear_believe61 fear_believe62 fear_believe63 fear_believe64 fear_believe65 fear_believe66 fear_believe67 fear_believe68 fear_believe69 fear_believe70 fear_believe71 fear_believe72 fear_believe73 fear_believe74 fear_believe75 fear_believe76 fear_believe77 fear_believe78 fear_believe79 fear_believe80 fear_believe81 fear_believe82 fear_believe83 fear_believe84 fear_believe85 fear_believe86 fear_believe87 fear_believe88 fear_believe89 fear_believe90 fear_believe91 fear_believe92 fear_believe93 fear_believe94 fear_believe95;




surprise_believe1 = EyeOpening(5) * MouthOpening(5) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1);
surprise_believe2 = EyeOpening(5) * MouthOpening(5) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe3 = EyeOpening(5) * MouthOpening(5) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1);
surprise_believe4 = EyeOpening(5) * MouthOpening(5) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe5 = EyeOpening(4) * MouthOpening(5) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1);
surprise_believe6 = EyeOpening(4) * MouthOpening(5) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe7 = EyeOpening(4) * MouthOpening(5) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1);
surprise_believe8 = EyeOpening(4) * MouthOpening(5) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe9 = EyeOpening(5) * MouthOpening(4) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe10 = EyeOpening(5) * MouthOpening(4) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe11 = EyeOpening(5) * MouthOpening(4) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe12 = EyeOpening(5) * MouthOpening(4) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe13 = EyeOpening(4) * MouthOpening(4) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe14 = EyeOpening(4) * MouthOpening(4) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe15 = EyeOpening(4) * MouthOpening(4) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe16 = EyeOpening(4) * MouthOpening(4) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe17 = EyeOpening(5) * MouthOpening(2) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe18 = EyeOpening(4) * MouthOpening(2) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(1) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe19 = EyeOpening(4) * MouthOpening(2) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe20 = EyeOpening(5) * MouthOpening(2) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(3) * ((1-NoseSideWrinkle)) * MouthLength(1) * MeanIntensity(2);
surprise_believe21 = EyeOpening(5) * MouthOpening(3) * (1-MouthCornerDisplacement(6)) * max(MouthCornerDisplacement) * EyebrowConstriction(2) * ((1-NoseSideWrinkle)) * MeanIntensity(2);

very_surprise = max([surprise_believe1 surprise_believe2 surprise_believe3 surprise_believe4 surprise_believe5 surprise_believe6 surprise_believe7 surprise_believe8  surprise_believe9 surprise_believe10 surprise_believe11 surprise_believe12 surprise_believe13 surprise_believe14 surprise_believe15 surprise_believe16 surprise_believe21]);
moderately_surprise = max([surprise_believe17 surprise_believe18 surprise_believe19 surprise_believe20]);

[surprise_believe surprise_intensity] = max([very_surprise moderately_surprise]);
% [surprise_believe sui] = max([surprise_believe1 surprise_believe2 surprise_believe3 surprise_believe4 surprise_believe5 surprise_believe6 surprise_believe7 surprise_believe8 surprise_believe9 surprise_believe10 surprise_believe11 surprise_believe12 surprise_believe13 surprise_believe14 surprise_believe15 surprise_believe16 surprise_believe17 surprise_believe18 surprise_believe19 surprise_believe20 surprise_believe21;



angry_believe1 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
angry_believe2 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
angry_believe3 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
angry_believe4 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
angry_believe5 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
angry_believe6 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
angry_believe7 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
angry_believe8 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * SadBelieve(1);
angry_believe9 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
angry_believe10 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * SadBelieve(1);
angry_believe11 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(4);
angry_believe12 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * SadBelieve(1);
angry_believe13 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
angry_believe14 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
angry_believe15 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
angry_believe16 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * SadBelieve(1);
angry_believe17 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
angry_believe18 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * SadBelieve(1);
angry_believe19 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
angry_believe20 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * SadBelieve(1);
angry_believe21 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
angry_believe22 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * SadBelieve(1);
angry_believe23 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
angry_believe24 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * SadBelieve(1);
% angry_believe25 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
% angry_believe26 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
% angry_believe27 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
% angry_believe28 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
% angry_believe29 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
% angry_believe30 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
% angry_believe31 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
% angry_believe32 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
% angry_believe33 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
% angry_believe34 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
% angry_believe35 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(4);
% angry_believe36 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(4);
% angry_believe37 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
% angry_believe38 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
% angry_believe39 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
% angry_believe40 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
% angry_believe41 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
% angry_believe42 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
% angry_believe43 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
% angry_believe44 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
% angry_believe45 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
% angry_believe46 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(4);
% angry_believe47 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
% angry_believe48 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
% angry_believe49 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
% angry_believe50 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
% angry_believe51 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
% angry_believe52 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
% angry_believe53 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
% angry_believe54 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
% angry_believe55 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
% angry_believe56 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
% angry_believe57 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
% angry_believe58 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
% angry_believe59 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
% angry_believe60 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
% angry_believe61 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
angry_believe62 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * SadBelieve(1);

very_angry = max([angry_believe1 angry_believe2 angry_believe3 angry_believe4 angry_believe5 angry_believe6 angry_believe7 angry_believe8 angry_believe9 angry_believe11 angry_believe13 angry_believe14 angry_believe15 angry_believe16 angry_believe17 angry_believe18 angry_believe19 angry_believe21 angry_believe23 angry_believe62]);
moderately_angry = max([angry_believe10 angry_believe12 angry_believe20 angry_believe22 angry_believe24]);

[angry_believe angry_intensity] = max([very_angry moderately_angry]);
% [angry_believe ai] = max([angry_believe1 angry_believe2 angry_believe3 angry_believe4 angry_believe5 angry_believe6 angry_believe7 angry_believe8 angry_believe9 angry_believe10 angry_believe11 angry_believe12 angry_believe13 angry_believe14 angry_believe15 angry_believe16 angry_believe17 angry_believe18 angry_believe19 angry_believe20 angry_believe21 angry_believe22 angry_believe23 angry_believe24 angry_believe25 angry_believe26 angry_believe27 angry_believe28 angry_believe29 angry_believe30 angry_believe31 angry_believe32 angry_believe33 angry_believe34 angry_believe35 angry_believe36 angry_believe37 angry_believe38 angry_believe39 angry_believe40 angry_believe41 angry_believe42 angry_believe43 angry_believe44 angry_believe45 angry_believe46 angry_believe47 angry_believe48 angry_believe49 angry_believe50 angry_believe51 angry_believe52 angry_believe53 angry_believe54 angry_believe55 angry_believe56 angry_believe57 angry_believe58 angry_believe59 angry_believe60 angry_believe61;
% [angry_believe ai] = max([angry_believe1 angry_believe2 angry_believe3 angry_believe4 angry_believe5 angry_believe6 angry_believe7 angry_believe8 angry_believe9 angry_believe10 angry_believe11 angry_believe12 angry_believe13 angry_believe14 angry_believe15 angry_believe16 angry_believe17 angry_believe18 angry_believe19 angry_believe20 angry_believe21 angry_believe22 angry_believe23 angry_believe24 angry_believe62;



disgust_believe1 = EyeOpening(1) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe2 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe3 = EyeOpening(1) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe4 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe5 = EyeOpening(1) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe6 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe7 = EyeOpening(1) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe8 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe9 = EyeOpening(1) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe10 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe11 = EyeOpening(1) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe12 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe13 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe14 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe15 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe16 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe17 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe18 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe19 = EyeOpening(1) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe20 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe21 = EyeOpening(1) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe22 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe23 = EyeOpening(1) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe24 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe25 = EyeOpening(1) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe26 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe27 = EyeOpening(1) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe28 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe29 = EyeOpening(1) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe30 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe31 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe32 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe33 = EyeOpening(3) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe34 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe35 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe36 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe37 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe38 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe39 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe40 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe41 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe42 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe43 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe44 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe45 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe46 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe47 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe48 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe49 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe50 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe51 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe52 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe53 = EyeOpening(1) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
disgust_believe54 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
disgust_believe55 = EyeOpening(1) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
disgust_believe56 = EyeOpening(2) * MouthOpening(2) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
disgust_believe57 = EyeOpening(1) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
disgust_believe58 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
disgust_believe59 = EyeOpening(1) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
disgust_believe60 = EyeOpening(2) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
disgust_believe61 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
disgust_believe62 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
disgust_believe63 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(4);
disgust_believe64 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(6) * EyebrowConstriction(5);
disgust_believe65 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe66 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe67 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe68 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe69 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe70 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe71 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe72 = EyeOpening(4) * MouthOpening(2) * MouthCornerDisplacement(5) * EyebrowConstriction(5);
disgust_believe73 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe74 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe75 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(5) * EyebrowConstriction(4);
disgust_believe76 = EyeOpening(4) * MouthOpening(3) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe77 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe78 = EyeOpening(1) * MouthOpening(3) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe79 = EyeOpening(1) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe80 = EyeOpening(1) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe81 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe82 = EyeOpening(2) * MouthOpening(3) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe83 = EyeOpening(2) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe84 = EyeOpening(2) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe85 = EyeOpening(1) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe86 = EyeOpening(1) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe87 = EyeOpening(2) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe88 = EyeOpening(2) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe89 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe90 = EyeOpening(3) * MouthOpening(3) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe91 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe92 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe93 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe94 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe95 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe96 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(5);
disgust_believe97 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe98 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe99 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * NoseSideWrinkle * MouthLength(1);
disgust_believe100 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe101 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle * MouthLength(1);
disgust_believe102 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe103 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe104 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe105 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe106 = EyeOpening(3) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe107 = EyeOpening(4) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * NoseSideWrinkle * MouthLength(1);
disgust_believe108 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe109 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe110 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe111 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe112 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe113 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe114 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe115 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe116 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe117 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe118 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe119 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe120 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe121 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe122 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe123 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe124 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe125 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * NoseSideWrinkle * MouthLength(1);
disgust_believe126 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe127 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe128 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe129 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe130 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe131 = EyeOpening(5) * MouthOpening(1) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe132 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe133 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe134 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe135 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe136 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe137 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe138 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe139 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe140 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe141 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe142 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe143 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe144 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe145 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe146 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(5) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe147 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe148 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe149 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe150 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(4) * NoseSideWrinkle;
disgust_believe151 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(6) * EyebrowConstriction(5) * NoseSideWrinkle;
disgust_believe152 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * NoseSideWrinkle * MouthLength(1);
disgust_believe153 = EyeOpening(1) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * NoseSideWrinkle * MouthLength(1);
disgust_believe154 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * NoseSideWrinkle * MouthLength(1);
disgust_believe155 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * NoseSideWrinkle * MouthLength(1);
disgust_believe156 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * NoseSideWrinkle * MouthLength(1);
disgust_believe157 = EyeOpening(3) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * NoseSideWrinkle * MouthLength(1);
disgust_believe158 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(4);
disgust_believe159 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(5);
disgust_believe160 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe161 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe162 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe163 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(3) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe164 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(3) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe165 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe166 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe167 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(3) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe168 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(3) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe169 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe170 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe171 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe172 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe173 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe174 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(3) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe175 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(3) * EyebrowConstriction(4) * MouthLength(1);
disgust_believe176 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe177 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(2) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe178 = EyeOpening(5) * MouthOpening(4) * MouthCornerDisplacement(3) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe179 = EyeOpening(5) * MouthOpening(5) * MouthCornerDisplacement(3) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe180 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(5) * MouthLength(1);
disgust_believe181 = EyeOpening(2) * MouthOpening(5) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * MouthLength(1);
disgust_believe182 = EyeOpening(3) * MouthOpening(4) * MouthCornerDisplacement(4) * EyebrowConstriction(2) * MouthLength(1);
disgust_believe183 = EyeOpening(4) * MouthOpening(4) * MouthCornerDisplacement(1) * EyebrowConstriction(2) * NoseSideWrinkle * MouthLength(1);
disgust_believe184 = EyeOpening(4) * MouthOpening(5) * MouthCornerDisplacement(4) * EyebrowConstriction(3) * NoseSideWrinkle * MouthLength(1);

very_disgust = max([disgust_believe1 disgust_believe2 disgust_believe3 disgust_believe4 disgust_believe5 disgust_believe6 disgust_believe7 disgust_believe8 disgust_believe9 disgust_believe10 disgust_believe11 disgust_believe12 disgust_believe13 disgust_believe14 disgust_believe15 disgust_believe16 disgust_believe17 disgust_believe18 disgust_believe19 disgust_believe20 disgust_believe21 disgust_believe22 disgust_believe23 disgust_believe24 disgust_believe25 disgust_believe26 disgust_believe27 disgust_believe28 disgust_believe29 disgust_believe30 disgust_believe31 disgust_believe32 disgust_believe33 disgust_believe34 disgust_believe35 disgust_believe36 disgust_believe37 disgust_believe38 disgust_believe39 disgust_believe40 disgust_believe41 disgust_believe42 disgust_believe43 disgust_believe46 disgust_believe47 disgust_believe48 disgust_believe49 disgust_believe50 disgust_believe51 disgust_believe52 disgust_believe53 disgust_believe54 disgust_believe54 disgust_believe55 disgust_believe56 disgust_believe57 disgust_believe58 disgust_believe58 disgust_believe59 disgust_believe60 disgust_believe61 disgust_believe62 disgust_believe63 disgust_believe64 disgust_believe68 disgust_believe69 disgust_believe70 disgust_believe71 disgust_believe72 disgust_believe73 disgust_believe74 disgust_believe75 disgust_believe76 disgust_believe77 disgust_believe78 disgust_believe79 disgust_believe80 disgust_believe81 disgust_believe82 disgust_believe83 disgust_believe84 disgust_believe85 disgust_believe86 disgust_believe87 disgust_believe88 disgust_believe89 disgust_believe90 disgust_believe91 disgust_believe92 disgust_believe93 disgust_believe94 disgust_believe95 disgust_believe96 disgust_believe97 disgust_believe98 disgust_believe99 disgust_believe100 disgust_believe101 disgust_believe102 disgust_believe103 disgust_believe104 disgust_believe105 disgust_believe106 disgust_believe107 disgust_believe108 disgust_believe109 disgust_believe110 disgust_believe111 disgust_believe112 disgust_believe113 disgust_believe114 disgust_believe115 disgust_believe116 disgust_believe117 disgust_believe118 disgust_believe119 disgust_believe120 disgust_believe121 disgust_believe122 disgust_believe123 disgust_believe124 disgust_believe125 disgust_believe126 disgust_believe127 disgust_believe128 disgust_believe129 disgust_believe130 disgust_believe131 disgust_believe132 disgust_believe133 disgust_believe134 disgust_believe135 disgust_believe136 disgust_believe137 disgust_believe138 disgust_believe139 disgust_believe140 disgust_believe141 disgust_believe142 disgust_believe143 disgust_believe144 disgust_believe145 disgust_believe146 disgust_believe147 disgust_believe148 disgust_believe149 disgust_believe150 disgust_believe151 disgust_believe152 disgust_believe153 disgust_believe154 disgust_believe155 disgust_believe156 disgust_believe157 disgust_believe158 disgust_believe159 disgust_believe160 disgust_believe161 disgust_believe162 disgust_believe163 disgust_believe164 disgust_believe165 disgust_believe166 disgust_believe167 disgust_believe168 disgust_believe169 disgust_believe170 disgust_believe171 disgust_believe172 disgust_believe172 disgust_believe173 disgust_believe174 disgust_believe175 disgust_believe176 disgust_believe177 disgust_believe178 disgust_believe179 disgust_believe180 disgust_believe181 disgust_believe182 disgust_believe183  disgust_believe184]);
moderately_disgust = max([disgust_believe44 disgust_believe45 disgust_believe65 disgust_believe67]);
not_so_disgust = disgust_believe66;

[disgust_believe disgust_intensity] = max([very_disgust moderately_disgust not_so_disgust]);
% [disgust_believe di] = max([disgust_believe1 disgust_believe2 disgust_believe3 disgust_believe4 disgust_believe5 disgust_believe6 disgust_believe7 disgust_believe8 disgust_believe9 disgust_believe10 disgust_believe11 disgust_believe12 disgust_believe13 disgust_believe14 disgust_believe15 disgust_believe16 disgust_believe17 disgust_believe18 disgust_believe19 disgust_believe20 disgust_believe21 disgust_believe22 disgust_believe23 disgust_believe24 disgust_believe25 disgust_believe26 disgust_believe27 disgust_believe28 disgust_believe29 disgust_believe30 disgust_believe31 disgust_believe32 disgust_believe33 disgust_believe34 disgust_believe35 disgust_believe36 disgust_believe37 disgust_believe38 disgust_believe39 disgust_believe40 disgust_believe41 disgust_believe42 disgust_believe43 disgust_believe44 disgust_believe45 disgust_believe46 disgust_believe47 disgust_believe48 disgust_believe49 disgust_believe50 disgust_believe51 disgust_believe52 disgust_believe53 disgust_believe54 disgust_believe54 disgust_believe55 disgust_believe56 disgust_believe57 disgust_believe58 disgust_believe58 disgust_believe59 disgust_believe60 disgust_believe61 disgust_believe62 disgust_believe63 disgust_believe64 disgust_believe65 disgust_believe66 disgust_believe67 disgust_believe68 disgust_believe69 disgust_believe70 disgust_believe71 disgust_believe72 disgust_believe73 disgust_believe74 disgust_believe75 disgust_believe76 disgust_believe77 disgust_believe78 disgust_believe79 disgust_believe80 disgust_believe81 disgust_believe82 disgust_believe83 disgust_believe84 disgust_believe85 disgust_believe86 disgust_believe87 disgust_believe88 disgust_believe89 disgust_believe90 disgust_believe91 disgust_believe92 disgust_believe93 disgust_believe94 disgust_believe95 disgust_believe96 disgust_believe97 disgust_believe98 disgust_believe99 disgust_believe100 disgust_believe101 disgust_believe102 disgust_believe103 disgust_believe104 disgust_believe105 disgust_believe106 disgust_believe107 disgust_believe108 disgust_believe109 disgust_believe110 disgust_believe111 disgust_believe112 disgust_believe113 disgust_believe114 disgust_believe115 disgust_believe116 disgust_believe117 disgust_believe118 disgust_believe119 disgust_believe120 disgust_believe121 disgust_believe122 disgust_believe123 disgust_believe124 disgust_believe125 disgust_believe126 disgust_believe127 disgust_believe128 disgust_believe129 disgust_believe130 disgust_believe131 disgust_believe132 disgust_believe133 disgust_believe134 disgust_believe135 disgust_believe136 disgust_believe137 disgust_believe138 disgust_believe139 disgust_believe140 disgust_believe141 disgust_believe142 disgust_believe143 disgust_believe144 disgust_believe145 disgust_believe146 disgust_believe147 disgust_believe148 disgust_believe149 disgust_believe150 disgust_believe151 disgust_believe152 disgust_believe153 disgust_believe154 disgust_believe155 disgust_believe156 disgust_believe157 disgust_believe158 disgust_believe159 disgust_believe160 disgust_believe161 disgust_believe162 disgust_believe163 disgust_believe164 disgust_believe165 disgust_believe166 disgust_believe167 disgust_believe168 disgust_believe169 disgust_believe170 disgust_believe171 disgust_believe172 disgust_believe172 disgust_believe173 disgust_believe174 disgust_believe175 disgust_believe176 disgust_believe177 disgust_believe178 disgust_believe179 disgust_believe180 disgust_believe181 disgust_believe182 disgust_believe183  disgust_believe184;



EmotionBelieve = [happy_believe sad_believe fear_believe surprise_believe angry_believe disgust_believe];
[Max_Believe max_index] = max(EmotionBelive);

if(length(Max_Believe) > 1)
    Emotion_Intencity = 'Blend Emotions';
else
    switch(max_index)
        case 1
            if(happy_intensity == 1)
                Emotion_Intencity = 'Very Happy';
            else
                if(happy_intensity == 2)
                    Emotion_Intencity = 'Moderately Happy';
                else
                    Emotion_Intencity = 'Not-So Happy';
                end
            end
            
        case 2
            if(sad_intensity == 1)
                Emotion_Intencity = 'Very Sad';
            else
                if(sad_intensity == 2)
                    Emotion_Intencity = 'Moderately Sad';
                else
                    Emotion_Intencity = 'Not-So Sad';
                end
            end

        case 3
            if(fear_intensity == 1)
                Emotion_Intencity = 'Very Fear';
            else
                if(fear_intensity == 2)
                    Emotion_Intencity = 'Moderately Fear';
                else
                    Emotion_Intencity = 'Not-So Fear';
                end
            end
            
        case 4
            if(surprise_intensity == 1)
                Emotion_Intencity = 'Very Surprise';
            else
                Emotion_Intencity = 'Moderately Surprise';
            end
            
        case 5
            if(angry_intensity == 1)
                Emotion_Intencity = 'Very Angry';
            else
                Emotion_Intencity = 'Moderately Angry';
            end
            
        case 6
            if(disgust_intensity == 1)
                Emotion_Intencity = 'Very Disgust';
            else
                if(disgust_intensity == 2)
                    Emotion_Intencity = 'Moderately Disgust';
                else
                    Emotion_Intencity = 'Not-So Disgust';
                end
            end

    end
end




function txtProcessTime_Callback(hObject, eventdata, handles)
% hObject    handle to txtProcessTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtProcessTime as text
%        str2double(get(hObject,'String')) returns contents of txtProcessTime as a double


% --- Executes during object creation, after setting all properties.
function txtProcessTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtProcessTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowFacialFeatures.
function ShowFacialFeatures_Callback(hObject, eventdata, handles)
% hObject    handle to ShowFacialFeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowFacialFeatures


