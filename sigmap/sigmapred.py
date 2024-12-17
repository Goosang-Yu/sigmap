import joblib, inspect
import pandas as pd
from time import time
from sklearn.preprocessing import MinMaxScaler
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# 새로 정리한 코드들
import sigmap
from sigmap.setup import readseq, scanner, designer
from sigmap.feature import FeatureExtraction

class SigmaFactor:
    def __init__(self, ):
        
        self.start = time()
        pass
    
    def predict(self, inFile:str, threshold=0.34):
        
        self.s70 = readseq(inFile)
        feat = FeatureExtraction(self.s70)
        
        self.list_seq   = feat.list_seq
        self.list_seqID = feat.list_seqID
        self.data       = feat.cdk
        
        result = self._run_model(self.data, threshold=threshold)
        
        return result
        
        
    def scan(self, inFile:str, threshold=0.34):
        
        self.s70 = scanner(inFile)
        feat = FeatureExtraction(self.s70)
        
        self.list_seq   = feat.list_seq
        self.list_seqID = feat.list_seqID
        self.data       = feat.cdk
        
        result = self._run_model(self.data, threshold=threshold)
        
        return result

    def design(self, inFile:str, threshold=0.34):
        
        self.s70 = designer(inFile)
        feat = FeatureExtraction(self.s70)
        
        self.list_seq   = feat.list_seq
        self.list_seqID = feat.list_seqID
        self.data       = feat.cdk
        
        result = self._run_model(self.data, threshold=threshold)
        
        return result

    def _run_model(self, final_df, threshold=0.34):

        # Load model
        base_dir  = inspect.getfile(sigmap).replace('__init__.py', '')
        model_dir = f'{base_dir}/model/svc_model_sigma70'
        model     = joblib.load(model_dir)
        
        # Features scaling
        to_scale = final_df.values
        scaler   = MinMaxScaler()
        scaler.fit(to_scale)
        scaled_data = scaler.transform(to_scale)

        # Prediction
        x_train = pd.DataFrame(scaled_data)
        y_pred  = model.predict_proba(x_train)
        y_class = []

        # Assigning class using threshold
        for i in y_pred:
            if(i[1] > threshold):
                y_class.append('Promoter')
            else:
                y_class.append('Non-Promoter')
                
        result = [[self.list_seqID[i],self.list_seq[i],round(val[1],3),y_class[i]] for i, val in enumerate(y_pred)]
        
        df_out = pd.DataFrame(result)
        df_out.columns = ["ID","Sequence"," Score", "Prediction"]
        
        return df_out
    
