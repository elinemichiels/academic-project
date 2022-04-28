import pickle
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble  import AdaBoostClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.preprocessing import StandardScaler

class QuinielaModel:

    def __init__(self):
        """Defines the classifiers that form part of the vote as well as the
        VotingClassifier object.
        """

        self.knn = KNeighborsClassifier(n_neighbors=8, weights='distance' )
        self.dt = DecisionTreeClassifier(max_depth=3, class_weight='balanced', random_state=20)
        self.dt_2 = DecisionTreeClassifier(max_depth=3, class_weight='balanced', random_state=10)
        self.adb = AdaBoostClassifier(base_estimator=self.dt_2 , n_estimators=400)
        self.classifiers = [('k nearest', self.knn),
              ('Decision Tree', self.dt),
              ('Ada', self.adb)]
        self.vc = VotingClassifier(estimators=self.classifiers)

        self.features = ['season','matchday','home_rank','away_rank','home_form','away_form','home_GF_pg',
              'home_GA_pg','away_GF_pg','away_GA_pg','home_rank_HT','away_rank_AT','last_conf']
        
        self.scaler = StandardScaler()

    def train(self, train_data):
        """Fits the VotingClassifier model after scaling the train_data
        :param train_data: DataFrame with the matches to be used for training
        """

        train_data = train_data[train_data['season']<2021]
        if train_data.empty:
            raise ValueError(f"Please specify at least a training season earlier than 2021-2022.")

        X = train_data[self.features].copy()
        y = train_data['result'].copy()

        self.scaler.fit(X)
        X_scaled = self.scaler.transform(X)
        self.vc.fit(X_scaled, y)

    def predict(self, predict_data):
        """Predicts the matches in predict_data using the self.vc model.
        :param predict_data: DataFrame with the matches to predict
        """
        
        X_pred = predict_data[self.features]
        X_pred_scaled = self.scaler.transform(X_pred)
        return self.vc.predict(X_pred_scaled)

    @classmethod
    def load(cls, filename):
        """ Load model from file """
        with open(filename, "rb") as f:
            model = pickle.load(f)
            assert type(model) == cls
        return model

    def save(self, filename):
        """ Save a model in a file """
        with open(filename, "wb") as f:
            pickle.dump(self, f)
