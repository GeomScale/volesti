function mpt_kblink(kb_id, varargin)
% Shows a link to a corresponding knowledge-base article

% TODO: keep in sync with the MPT3 website
base_url = 'http://control.ee.ethz.ch/~mpt/redmine/knowledgebase/articles/';
kb_url = [base_url num2str(kb_id)];
fprintf('\nSee <a href="%s">the corresponding knowledge base article</a> for the explanation of this error.\n', kb_url);

end
